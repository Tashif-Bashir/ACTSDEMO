#!/usr/bin/env python3
import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    addCKFTracks,
    CKFPerformanceConfig,
    addVertexFitting,
    VertexFinder,
)

ttbar_pu200 = False
u = acts.UnitConstants
geo_dir = pathlib.Path("acts-itk")
outputDir = pathlib.Path.cwd() / "itk_output"
outputDir.mkdir(exist_ok=True)
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector, trackingGeometry, decorators = acts.examples.itk.buildITkGeometry(geo_dir)
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=100, numThreads=-1)

if not ttbar_pu200:
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, uniform=True),
        ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )
else:
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )

addFatras(
    s,
    trackingGeometry,
    field,
    ParticleSelectorConfig(eta=(-4.0, 4.0), pt=(150 * u.MeV, None), removeNeutral=True)
    if ttbar_pu200
    else ParticleSelectorConfig(),
    outputDirRoot=outputDir,
    rnd=rnd,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir / "itk-hgtd/itk-smearing-config.json",
    outputDirRoot=outputDir,
    rnd=rnd,
)

addSeeding(
    s,
    trackingGeometry,
    field,
    TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-4.0, 4.0), nHits=(9, None))
    if ttbar_pu200
    else TruthSeedRanges(),
    *acts.examples.itk.itkSeedingAlgConfig("PixelSpacePoints"),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json",
    outputDirRoot=outputDir,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    CKFPerformanceConfig(ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=6),
    outputDirRoot=outputDir,
)

s.addAlgorithm(
    acts.examples.TrackSelector(
        level=acts.logging.INFO,
        inputTrackParameters="fittedTrackParameters",
        outputTrackParameters="trackparameters",
        outputTrackIndices="outputTrackIndices",
        removeNeutral=True,
        absEtaMax=4.0,
        loc0Max=4.0 * u.mm,  # rho max
        ptMin=1.0 * u.GeV,
    )
)

addVertexFitting(
    s,
    field,
    vertexFinder=VertexFinder.Iterative,
    outputDirRoot=outputDir,
    trajectories="trajectories",
)

s.run()
