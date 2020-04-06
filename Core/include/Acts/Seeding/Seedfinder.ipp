// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <cmath>
#include <numeric>
#include <type_traits>
#include <algorithm>
#include <chrono>
#include <numeric>

namespace Acts {

  template <typename external_spacepoint_t, typename platform_t>
  Seedfinder<external_spacepoint_t, platform_t>::Seedfinder(
    Acts::SeedfinderConfig<external_spacepoint_t> config)
    : m_config(std::move(config)) {
  // calculation of scattering using the highland formula
  // convert pT to p once theta angle is known
  m_config.highland = 13.6 * std::sqrt(m_config.radLengthPerSeed) *
                      (1 + 0.038 * std::log(m_config.radLengthPerSeed));
  float maxScatteringAngle = m_config.highland / m_config.minPt;
  m_config.maxScatteringAngle2 = maxScatteringAngle * maxScatteringAngle;

  // helix radius in homogeneous magnetic field. Units are Kilotesla, MeV and
  // millimeter
  // TODO: change using ACTS units
  m_config.pTPerHelixRadius = 300. * m_config.bFieldInZ;
  m_config.minHelixDiameter2 =
      std::pow(m_config.minPt * 2 / m_config.pTPerHelixRadius, 2);
  m_config.pT2perRadius =
      std::pow(m_config.highland / m_config.pTPerHelixRadius, 2);    
  }

  
  template< typename external_spacepoint_t, typename platform_t>
  template< typename T, typename sp_range_t>
  typename std::enable_if< std::is_same<T, Acts::CPU>::value, std::vector<Seed<external_spacepoint_t> > >::type
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  int i_m=0;
  for (auto spM : middleSPs) {    
    i_m++;
    
    float rM = spM->radius();
    float zM = spM->z();
    float varianceRM = spM->varianceR();
    float varianceZM = spM->varianceZ();

    // Doublet search    
    auto compatBottomSP =
      SeedfinderCPUFunctions<external_spacepoint_t,
			     sp_range_t>::searchDoublet(true, bottomSPs, *spM, m_config);
    
    // no bottom SP found -> try next spM
    if (compatBottomSP.empty()) {
      continue;
    }

    auto compatTopSP =
      SeedfinderCPUFunctions<external_spacepoint_t,
			     sp_range_t>::searchDoublet(false, topSPs, *spM, m_config);

    // no top SP found -> try next spM
    if (compatTopSP.empty()) {
      continue;
    }

    // contains parameters required to calculate circle with linear equation
    
    // ...for bottom-middle
    std::vector<LinCircle> linCircleBottom;
    // ...for middle-top
    std::vector<LinCircle> linCircleTop;
    
    SeedfinderCPUFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatBottomSP, *spM, true, linCircleBottom);
    SeedfinderCPUFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatTopSP, *spM, false, linCircleTop);

    auto seedsPerSpM = SeedfinderCPUFunctions<external_spacepoint_t,sp_range_t>::searchTriplet(*spM, compatBottomSP, compatTopSP, linCircleBottom, linCircleTop, m_config);
    m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);   
  }
  
  return outputVec;
  }

#ifdef  ACTS_HAS_CUDA
  
  // CUDA seed finding
  template< typename external_spacepoint_t, typename platform_t>
  template< typename T, typename sp_range_t>
  typename std::enable_if< std::is_same<T, Acts::CUDA>::value, std::vector<Seed<external_spacepoint_t> > >::type
  Seedfinder<external_spacepoint_t, platform_t>::createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const {
  std::vector<Seed<external_spacepoint_t>> outputVec;

  CUDAArray<unsigned char> isBottom_cuda(1);

  unsigned char true_cpu  = true;
  unsigned char false_cpu = false;
  CUDAArray<unsigned char> true_cuda(1,&true_cpu,1);  
  CUDAArray<unsigned char> false_cuda(1,&false_cpu,1);    
  CUDAArray<float> deltaRMin_cuda(1, &m_config.deltaRMin, 1);
  CUDAArray<float> deltaRMax_cuda(1, &m_config.deltaRMax, 1);
  CUDAArray<float> cotThetaMax_cuda(1, &m_config.cotThetaMax, 1);
  CUDAArray<float> collisionRegionMin_cuda(1, &m_config.collisionRegionMin, 1);
  CUDAArray<float> collisionRegionMax_cuda(1, &m_config.collisionRegionMax, 1);  
  CUDAArray<float> maxScatteringAngle2_cuda(1, &m_config.maxScatteringAngle2,1);
  CUDAArray<float> sigmaScattering_cuda(1, &m_config.sigmaScattering,1);
  CUDAArray<float> minHelixDiameter2_cuda(1, &m_config.minHelixDiameter2,1);
  CUDAArray<float> pT2perRadius_cuda(1, &m_config.pT2perRadius,1);
  CUDAArray<float> impactMax_cuda(1, &m_config.impactMax,1);
  
  /*----------------------------------
     Algorithm 0. Matrix Flattening 
  ----------------------------------*/

  // Get Size of spacepoints
  int nMiddle = 0;
  int nBottom = 0;
  int nTop    = 0;

  for (auto sp: middleSPs) nMiddle++;
  for (auto sp: bottomSPs) nBottom++;
  for (auto sp: topSPs)    nTop++;

  if (nMiddle == 0 || nBottom == 0 || nTop == 0) return outputVec;

  // Define Matrix and Do flattening
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > middleSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > bottomSPvec;
  std::vector< const Acts::InternalSpacePoint<external_spacepoint_t>* > topSPvec;
  
  CPUMatrix<float> spMmat_cpu(nMiddle, 6); // x y z r varR varZ
  CPUMatrix<float> spBmat_cpu(nBottom, 6);
  CPUMatrix<float> spTmat_cpu(nTop   , 6);
  
  size_t mIdx=0;
  for (auto sp: middleSPs){
    spMmat_cpu.SetEl(mIdx,0,sp->x());
    spMmat_cpu.SetEl(mIdx,1,sp->y());
    spMmat_cpu.SetEl(mIdx,2,sp->z());
    spMmat_cpu.SetEl(mIdx,3,sp->radius());
    spMmat_cpu.SetEl(mIdx,4,sp->varianceR());
    spMmat_cpu.SetEl(mIdx,5,sp->varianceZ());
    middleSPvec.push_back(sp);
    mIdx++;
  }
  CUDAMatrix<float> spMmat_cuda(nMiddle, 6, &spMmat_cpu);
  CUDAArray<int>    nSpM_cuda(1,&nMiddle,1);
  
  size_t bIdx=0;
  for (auto sp: bottomSPs){
    spBmat_cpu.SetEl(bIdx,0,sp->x());
    spBmat_cpu.SetEl(bIdx,1,sp->y());
    spBmat_cpu.SetEl(bIdx,2,sp->z());
    spBmat_cpu.SetEl(bIdx,3,sp->radius());
    spBmat_cpu.SetEl(bIdx,4,sp->varianceR());
    spBmat_cpu.SetEl(bIdx,5,sp->varianceZ());
    bottomSPvec.push_back(sp);
    bIdx++;
  }

  size_t tIdx=0;
  for (auto sp: topSPs){
    spTmat_cpu.SetEl(tIdx,0,sp->x());
    spTmat_cpu.SetEl(tIdx,1,sp->y());
    spTmat_cpu.SetEl(tIdx,2,sp->z());
    spTmat_cpu.SetEl(tIdx,3,sp->radius());
    spTmat_cpu.SetEl(tIdx,4,sp->varianceR());
    spTmat_cpu.SetEl(tIdx,5,sp->varianceZ());
    topSPvec.push_back(sp);
    tIdx++;    
  }

  /*------------------------------------
     Algorithm 1. Doublet Search (DS)
  ------------------------------------*/
  
  int  BlockSize;
  dim3 DS_BlockSize;
  dim3 DS_GridSize(nMiddle,1,1);

  CUDAArray<float> rM_cuda(nMiddle, spMmat_cpu.GetEl(0,3), nMiddle);
  CUDAArray<float> zM_cuda(nMiddle, spMmat_cpu.GetEl(0,2), nMiddle);
  CUDAArray<float> rB_cuda(nBottom, spBmat_cpu.GetEl(0,3), nBottom);    
  CUDAArray<float> zB_cuda(nBottom, spBmat_cpu.GetEl(0,2), nBottom);
  CUDAArray<int>   nBottom_cuda(1, &nBottom, 1);
  CUDAArray<float> rT_cuda(nTop,    spTmat_cpu.GetEl(0,3), nTop);    
  CUDAArray<float> zT_cuda(nTop,    spTmat_cpu.GetEl(0,2), nTop);  
  CUDAArray<int>   nTop_cuda(1, &nTop, 1);
  
  ///// For bottom space points
  CUDAMatrix<unsigned char> isCompatBottomMat_cuda(nBottom, nMiddle);

  int  offsetDS;
  offsetDS=0;
  while(offsetDS<nBottom){
    DS_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nBottom-offsetDS), 1,1);
    SeedfinderCUDAKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  true_cuda.Get(),
					  rM_cuda.Get(), zM_cuda.Get(),
					  nBottom_cuda.Get(), rB_cuda.Get(offsetDS), zB_cuda.Get(offsetDS), 
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatBottomMat_cuda.GetEl(offsetDS,0));
    offsetDS+=DS_BlockSize.x;
  }
  CPUMatrix<unsigned char>  isCompatBottomMat_cpu(nBottom, nMiddle, &isCompatBottomMat_cuda);
  
  ///// For top space points
  CUDAMatrix<unsigned char> isCompatTopMat_cuda(nTop, nMiddle);
  offsetDS=0;
  while(offsetDS<nTop){
    DS_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nTop-offsetDS), 1,1);    
    SeedfinderCUDAKernels::searchDoublet( DS_GridSize, DS_BlockSize,
					  false_cuda.Get(),
					  rM_cuda.Get(), zM_cuda.Get(),
					  nTop_cuda.Get(), rT_cuda.Get(offsetDS), zT_cuda.Get(offsetDS), 
					  deltaRMin_cuda.Get(), deltaRMax_cuda.Get(), 
					  cotThetaMax_cuda.Get(),
					  collisionRegionMin_cuda.Get(),collisionRegionMax_cuda.Get(),
					  isCompatTopMat_cuda.GetEl(offsetDS,0));
    offsetDS+=DS_BlockSize.x;
  }
  CPUMatrix<unsigned char>  isCompatTopMat_cpu(nTop, nMiddle, &isCompatTopMat_cuda);

  /* -----------------------------------------
     Algorithm 2. Transform Coordinates (TC)
  -------------------------------------------*/
  
  std::vector< std::tuple<int, std::vector< int >, std::vector< int > > > mCompIndex;

  CPUArray<int>  nBcompMax_cpu(1); nBcompMax_cpu[0] = 0;
  CPUArray<int>  nTcompMax_cpu(1); nTcompMax_cpu[0] = 0;

  for (int i_m=0; i_m<nMiddle; i_m++){
    std::vector< int > bIndex;
    for (int i=0; i<nBottom; i++){
      if (*isCompatBottomMat_cpu.GetEl(i,i_m)) bIndex.push_back(i);
    }
    if (bIndex.empty()) continue;    
    std::vector< int > tIndex;
    for (int i=0; i<nTop; i++){
      if (*isCompatTopMat_cpu.GetEl(i,i_m)) tIndex.push_back(i);
    }
    if (tIndex.empty()) continue;
    
    auto tup = std::make_tuple(i_m, bIndex, tIndex);
    mCompIndex.push_back(tup);
    nBcompMax_cpu[0] = fmax(bIndex.size(), nBcompMax_cpu[0]);
    nTcompMax_cpu[0] = fmax(tIndex.size(), nTcompMax_cpu[0]);
  }
  
  CPUMatrix<float>  spMcompMat_cpu(mCompIndex.size(), 6);
  CUDAMatrix<float> spMcompMat_cuda(mCompIndex.size(),6);
  CUDAArray<int>    nBcompMax_cuda(1,nBcompMax_cpu.Get(),1);
  CUDAArray<int>    nTcompMax_cuda(1,nTcompMax_cpu.Get(),1);

  CPUMatrix<float>  spBcompMat_cpu   (nBcompMax_cpu[0], mCompIndex.size()*6);
  CUDAMatrix<float> spBcompMat_cuda  (nBcompMax_cpu[0], mCompIndex.size()*6);
  CUDAMatrix<float> circBcompMat_cuda(nBcompMax_cpu[0], mCompIndex.size()*6);  
  CPUMatrix<float>  spTcompMat_cpu   (nTcompMax_cpu[0], mCompIndex.size()*6);  
  CUDAMatrix<float> spTcompMat_cuda  (nTcompMax_cpu[0], mCompIndex.size()*6);   
  CUDAMatrix<float> circTcompMat_cuda(nTcompMax_cpu[0], mCompIndex.size()*6);
  
  for (int i_m=0; i_m<spMcompMat_cpu.GetNRows(); i_m++){
    auto mIndex = std::get<0>(mCompIndex[i_m]);
    auto bIndex = std::get<1>(mCompIndex[i_m]);
    auto tIndex = std::get<2>(mCompIndex[i_m]);

    for (int i=0; i<6; i++){
      spMcompMat_cpu.SetEl(i_m,i,std::move(*spMmat_cpu.GetEl(mIndex,i)));
    }      
    //spMcompMat_cpu.SetEl(i_m,0,*spMmat_cpu.GetEl(mIndex,0));
    //spMcompMat_cpu.SetEl(i_m,1,*spMmat_cpu.GetEl(mIndex,1));
    //spMcompMat_cpu.SetEl(i_m,2,*spMmat_cpu.GetEl(mIndex,2));
    //spMcompMat_cpu.SetEl(i_m,3,*spMmat_cpu.GetEl(mIndex,3));
    //spMcompMat_cpu.SetEl(i_m,4,*spMmat_cpu.GetEl(mIndex,4));
    //spMcompMat_cpu.SetEl(i_m,5,*spMmat_cpu.GetEl(mIndex,5));
    
    for (int i_b=0; i_b<bIndex.size(); i_b++){
      for (int i=0; i<6; i++){
	spBcompMat_cpu.SetEl(i_b, i_m*6+i, std::move(*spBmat_cpu.GetEl(bIndex[i_b],i)));
      }            
      //spBcompMat_cpu.SetEl(i_b, i_m*6+0, *spBmat_cpu.GetEl(bIndex[i_b],0));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+1, *spBmat_cpu.GetEl(bIndex[i_b],1));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+2, *spBmat_cpu.GetEl(bIndex[i_b],2));	    
      //spBcompMat_cpu.SetEl(i_b, i_m*6+3, *spBmat_cpu.GetEl(bIndex[i_b],3));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+4, *spBmat_cpu.GetEl(bIndex[i_b],4));
      //spBcompMat_cpu.SetEl(i_b, i_m*6+5, *spBmat_cpu.GetEl(bIndex[i_b],5));
    }

    for (int i_t=0; i_t<tIndex.size(); i_t++){
      for (int i=0; i<6; i++){
	spTcompMat_cpu.SetEl(i_t, i_m*6+i, std::move(*spTmat_cpu.GetEl(tIndex[i_t],i)));
      }                  
      //spTcompMat_cpu.SetEl(i_t, i_m*6+0, *spTmat_cpu.GetEl(tIndex[i_t],0));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+1, *spTmat_cpu.GetEl(tIndex[i_t],1));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+2, *spTmat_cpu.GetEl(tIndex[i_t],2));	    
      //spTcompMat_cpu.SetEl(i_t, i_m*6+3, *spTmat_cpu.GetEl(tIndex[i_t],3));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+4, *spTmat_cpu.GetEl(tIndex[i_t],4));
      //spTcompMat_cpu.SetEl(i_t, i_m*6+5, *spTmat_cpu.GetEl(tIndex[i_t],5));
    }    
  }

  spMcompMat_cuda.CopyH2D(spMcompMat_cpu.GetEl(), spMcompMat_cpu.GetSize());
  spBcompMat_cuda.CopyH2D(spBcompMat_cpu.GetEl(), spBcompMat_cpu.GetSize());
  spTcompMat_cuda.CopyH2D(spTcompMat_cpu.GetEl(), spTcompMat_cpu.GetSize());

  dim3 TC_GridSize(spMcompMat_cpu.GetNRows(),1,1);
  dim3 TC_BlockSize;
  
  // For bottom-middle
  int offsetTC;
  offsetTC=0;
  while(offsetTC<nBcompMax_cpu[0]){
    TC_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nBcompMax_cpu[0]-offsetTC),1,1);    
    SeedfinderCUDAKernels::transformCoordinates(TC_GridSize, TC_BlockSize,
						true_cuda.Get(),
						spMcompMat_cuda.GetEl(0,0),
						nBcompMax_cuda.Get(),
						spBcompMat_cuda.GetEl(offsetTC,0),
						circBcompMat_cuda.GetEl(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
  
  // For middle-top 
  offsetTC=0;
  while(offsetTC<nTcompMax_cpu[0]){
    TC_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, nTcompMax_cpu[0]-offsetTC),1,1);
    SeedfinderCUDAKernels::transformCoordinates(TC_GridSize, TC_BlockSize,
						false_cuda.Get(),
						spMcompMat_cuda.GetEl(0,0),
						nTcompMax_cuda.Get(),
						spTcompMat_cuda.GetEl(offsetTC,0),
						circTcompMat_cuda.GetEl(offsetTC,0));    
    offsetTC+=TC_BlockSize.x;
  }
  
  int nMcomp = mCompIndex.size();  
  CUDAArray<int>  nMcomp_cuda(1,&nMcomp,1);
  
  // retreive middle-bottom doublet circ information
  CPUMatrix<float> circBcompMat_cpu(nBcompMax_cpu[0], mCompIndex.size()*6);
  circBcompMat_cpu.CopyD2H(circBcompMat_cuda.GetEl(0,0), circBcompMat_cuda.GetSize());  

  std::vector<int> offsetVec(m_config.offsetVecSize);
  std::iota (std::begin(offsetVec), std::end(offsetVec), 0); // Fill with 0, 1, ..., 99.
  for (auto& el: offsetVec) el = el*MAX_BLOCK_SIZE;
  CUDAArray<int> offsetVec_cuda(offsetVec.size(),&offsetVec[0],offsetVec.size());
    
  const int         nTopPassLimit = m_config.nTopPassLimit;
  CUDAArray<int>    nTopPassLimit_cuda(1, &nTopPassLimit, 1);
  
  CUDAMatrix<int>   nTopPass_cuda(nBcompMax_cpu[0], nMcomp);
  CUDAMatrix<int>   tPassIndex_cuda(nTopPassLimit, nBcompMax_cpu[0]*nMcomp); 
  CUDAMatrix<float> curvatures_cuda(nTopPassLimit, nBcompMax_cpu[0]*nMcomp);       
  CUDAMatrix<float> impactparameters_cuda(nTopPassLimit, nBcompMax_cpu[0]*nMcomp);

  CPUMatrix<int>    nTopPass_cpu(nBcompMax_cpu[0], nMcomp);
  CPUMatrix<int>    tPassIndex_cpu(nTopPassLimit, nBcompMax_cpu[0]*nMcomp); 
  CPUMatrix<float>  curvatures_cpu(nTopPassLimit, nBcompMax_cpu[0]*nMcomp);
  CPUMatrix<float>  impactparameters_cpu(nTopPassLimit, nBcompMax_cpu[0]*nMcomp);

  cudaStream_t cuStream;
  cudaStreamCreate(&cuStream);
   
  for (int i_m=0; i_m<=mCompIndex.size(); i_m++){
    
    //cudaDeviceSynchronize();
    cudaStreamSynchronize(cuStream);

    // -----------------------------------
    //  Algorithm 3. Triplet Search (TS)
    //------------------------------------
    
    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;
    
    std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;
    
    if (i_m < mCompIndex.size()){
          
      // SEED Filtering is done ASYNCHRONOUSLY against Triplet search
      // Need to call it again after last iteration
            
      // For triplets collected at previous iteration
      
      auto mIndex = std::get<0>(mCompIndex[i_m]);
      auto bIndex = std::get<1>(mCompIndex[i_m]);
      auto tIndex = std::get<2>(mCompIndex[i_m]);
      
      dim3 TS_GridSize(bIndex.size(),1,1);
      dim3 TS_BlockSize;
      
      int i_ts = 0;    
      while ( offsetVec[i_ts] < tIndex.size() ){
	
	TS_BlockSize = dim3(fmin(MAX_BLOCK_SIZE, tIndex.size()-offsetVec[i_ts] ),
			    1,1);
      	
	SeedfinderCUDAKernels::searchTriplet(TS_GridSize, TS_BlockSize,
					     offsetVec_cuda.Get(i_ts),
					     nMcomp_cuda.Get(),
					     spMcompMat_cuda.GetEl(i_m,0),
					     nBcompMax_cuda.Get(),
					     spBcompMat_cuda.GetEl(0,6*i_m),
					     nTcompMax_cuda.Get(),
					     spTcompMat_cuda.GetEl(offsetVec[i_ts],6*i_m),
					     circBcompMat_cuda.GetEl(0,6*i_m),
					     circTcompMat_cuda.GetEl(offsetVec[i_ts],6*i_m),
					     // Seed finder config
					     maxScatteringAngle2_cuda.Get(),
					     sigmaScattering_cuda.Get(),
					     minHelixDiameter2_cuda.Get(),
					     pT2perRadius_cuda.Get(),
					     impactMax_cuda.Get(),
					     nTopPassLimit_cuda.Get(),
					     // output
					     nTopPass_cuda.GetEl(0,i_m),
					     tPassIndex_cuda.GetEl(0,nBcompMax_cpu[0]*i_m),
					     curvatures_cuda.GetEl(0,nBcompMax_cpu[0]*i_m),
					     impactparameters_cuda.GetEl(0,nBcompMax_cpu[0]*i_m),
					     &cuStream
					     );
	i_ts++;
      }
      
      nTopPass_cpu.CopyD2H(nTopPass_cuda.GetEl(0,i_m),
			   nBcompMax_cpu[0],
			   nBcompMax_cpu[0]*i_m,
			   &cuStream);
      tPassIndex_cpu.CopyD2H(tPassIndex_cuda.GetEl(0,nBcompMax_cpu[0]*i_m),
			     nTopPassLimit*nBcompMax_cpu[0],
			     nTopPassLimit*nBcompMax_cpu[0]*i_m,
			     &cuStream);
      curvatures_cpu.CopyD2H(curvatures_cuda.GetEl(0,nBcompMax_cpu[0]*i_m),
			     nTopPassLimit*nBcompMax_cpu[0],
			     nTopPassLimit*nBcompMax_cpu[0]*i_m,
			     &cuStream); 
      impactparameters_cpu.CopyD2H(impactparameters_cuda.GetEl(0,nBcompMax_cpu[0]*i_m),
				   nTopPassLimit*nBcompMax_cpu[0],
				   nTopPassLimit*nBcompMax_cpu[0]*i_m,
				   &cuStream);
    }

    // --------------------------------
    //  Algorithm 4. Seed Filter (SF)
    // --------------------------------
    
    if (i_m > 0){
      seedsPerSpM.clear();
      auto middleIdx     = std::get<0>(mCompIndex[i_m-1]);
      auto compBottomIdx = std::get<1>(mCompIndex[i_m-1]);
      auto compTopIdx    = std::get<2>(mCompIndex[i_m-1]);
      
      for (int i_b=0; i_b<compBottomIdx.size(); i_b++){
	int nTpass = *(nTopPass_cpu.GetEl(i_b,i_m-1));

	if (nTpass==0) continue;	
	
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();      
	float Zob = *(circBcompMat_cpu.GetEl(i_b,(i_m-1)*6));

	std::vector< std::tuple< int, int, int > > indexVec;
	for(int i_t=0; i_t<nTpass; i_t++){
	  int g_tIndex = compTopIdx[*tPassIndex_cpu.GetEl(i_t,i_b+(i_m-1)*nBcompMax_cpu[0])];
	  indexVec.push_back(std::make_tuple(g_tIndex,i_t,i_b));
	}
	sort(indexVec.begin(), indexVec.end()); 
	
	for(auto el: indexVec){
	  auto g_tIndex = std::get<0>(el);
	  auto tId      = std::get<1>(el);
	  auto bId      = std::get<2>(el);
	  
	  tVec.push_back(topSPvec[g_tIndex]);
	  curvatures.push_back(*curvatures_cpu.GetEl(tId,bId+(i_m-1)*nBcompMax_cpu[0]));
	  impactParameters.push_back(*impactparameters_cpu.GetEl(tId,bId+(i_m-1)*nBcompMax_cpu[0]));
	}
	
	std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> sameTrackSeeds;
	sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(*bottomSPvec[compBottomIdx[i_b]],
									     *middleSPvec[middleIdx],
									     tVec,
									     curvatures,
									     impactParameters,Zob)); 
	seedsPerSpM.insert(seedsPerSpM.end(),
			   std::make_move_iterator(sameTrackSeeds.begin()),
			   std::make_move_iterator(sameTrackSeeds.end()));	      
      }
      m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);      
    }    
    
    /*
    std::vector<std::pair<
      float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> seedsPerSpM;
    
    // SEED Filtering is done ASYNCHRONOUSLY against Triplet search
    // Need to call it again after last iteration
    
    std::vector<const InternalSpacePoint<external_spacepoint_t> *> tVec;
    std::vector<float> curvatures;
    std::vector<float> impactParameters;
    
    // For triplets collected at previous iteration
    if (i_m > 0){
      seedsPerSpM.clear();
      auto middleIdx     = std::get<0>(mCompIndex[i_m-1]);
      auto compBottomIdx = std::get<1>(mCompIndex[i_m-1]);
      auto compTopIdx    = std::get<2>(mCompIndex[i_m-1]);
      
      for (int i_b=0; i_b<compBottomIdx.size(); i_b++){
	if (nTopPass_cpu[i_b]==0) continue;
	
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();      
	float Zob = *(circBcompMat_cpu.GetEl(i_b,i_m*6));

	std::vector< std::tuple< int, int, int > > indexVec;
	for(int i_t=0; i_t<nTopPass_cpu[i_b]; i_t++){
	  int g_tIndex = compTopIdx[*tPassIndex_cpu.GetEl(i_t,i_b)];
	  indexVec.push_back(std::make_tuple(g_tIndex,i_t,i_b));
	}
	sort(indexVec.begin(), indexVec.end()); 
	
	for(auto el: indexVec){
	  auto g_tIndex = std::get<0>(el);
	  auto tId      = std::get<1>(el);
	  auto bId      = std::get<2>(el);
	  
	  tVec.push_back(topSPvec[g_tIndex]);
	  curvatures.push_back(*curvatures_cpu.GetEl(tId,bId));
	  impactParameters.push_back(*impactparameters_cpu.GetEl(tId,bId));
	}
	
	std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> sameTrackSeeds;
	sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(*bottomSPvec[compBottomIdx[i_b]],
									     *middleSPvec[middleIdx],
									     tVec,
									     curvatures,
									     impactParameters,Zob)); 
	seedsPerSpM.insert(seedsPerSpM.end(),
			   std::make_move_iterator(sameTrackSeeds.begin()),
			   std::make_move_iterator(sameTrackSeeds.end()));	      
      }
      m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);      
    }
    
    nTopPass_cpu.CopyD2H(nTopPass_cuda.Get(), nBcompMax_cpu[0]);
    tPassIndex_cpu.CopyD2H(tPassIndex_cuda.GetEl(0,0),             nTopPassLimit*nBcompMax_cpu[0]);
    curvatures_cpu.CopyD2H(curvatures_cuda.GetEl(0,0),             nTopPassLimit*nBcompMax_cpu[0]);
    impactparameters_cpu.CopyD2H(impactparameters_cuda.GetEl(0,0), nTopPassLimit*nBcompMax_cpu[0]);
        
    // For the last iteration
    if (i_m == mCompIndex.size()-1 ){
      seedsPerSpM.clear();
      auto middleIdx     = std::get<0>(mCompIndex[i_m]);
      auto compBottomIdx = std::get<1>(mCompIndex[i_m]);
      auto compTopIdx    = std::get<2>(mCompIndex[i_m]);
      
      for (int i_b=0; i_b<compBottomIdx.size(); i_b++){
	if (nTopPass_cpu[i_b]==0) continue;
	
	tVec.clear();
	curvatures.clear();
	impactParameters.clear();      
	float Zob = *(circBcompMat_cpu.GetEl(i_b,i_m*6));

	std::vector< std::tuple< int, int, int > > indexVec;
	for(int i_t=0; i_t<nTopPass_cpu[i_b]; i_t++){
	  int g_tIndex = compTopIdx[*tPassIndex_cpu.GetEl(i_t,i_b)];
	  indexVec.push_back(std::make_tuple(g_tIndex,i_t,i_b));
	}
	sort(indexVec.begin(), indexVec.end()); 
	
	for(auto el: indexVec){
	  auto g_tIndex = std::get<0>(el);
	  auto tId      = std::get<1>(el);
	  auto bId      = std::get<2>(el);
	  
	  tVec.push_back(topSPvec[g_tIndex]);
	  curvatures.push_back(*curvatures_cpu.GetEl(tId,bId));
	  impactParameters.push_back(*impactparameters_cpu.GetEl(tId,bId));
	}
	
	std::vector<std::pair<float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>> sameTrackSeeds;
	sameTrackSeeds = std::move(m_config.seedFilter->filterSeeds_2SpFixed(*bottomSPvec[compBottomIdx[i_b]],
									     *middleSPvec[middleIdx],
									     tVec,
									     curvatures,
									     impactParameters,Zob)); 
	seedsPerSpM.insert(seedsPerSpM.end(),
			   std::make_move_iterator(sameTrackSeeds.begin()),
			   std::make_move_iterator(sameTrackSeeds.end()));	      
      }
      m_config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);            
    }    
    */
  }  
  
  return outputVec;  
  }

#endif
  
}// namespace Acts
