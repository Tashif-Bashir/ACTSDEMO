#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union
import enum

from acts.examples import (
    Sequencer,
    ParticleSelector,
    ParticleSmearing,
    RootParticleReader,
    RootTrajectorySummaryReader,
    TrackSelector,
    GenericDetector,
)

import acts

from acts import UnitConstants as u

import pythia8


from acts.examples.reconstruction import VertexFinder, addVertexFitting

def runVertexFitting(
    field,
    outputDir: Path,
    outputRoot: bool = True,
    inputParticlePath: Optional[Path] = None,
    inputTrackSummary: Path = None,
    vertexFinder: VertexFinder = VertexFinder.Truth,
    s=None,
):
    s = s or Sequencer(events=100, numThreads=-1)

    logger = acts.logging.getLogger("VertexFittingExample")

    rnd = acts.examples.RandomNumbers(seed=42)

    inputParticles = "particles_input"
    if inputParticlePath is None:
        logger.info("Generating particles using Pythia8")
        pythia8.addPythia8(s, rnd)
    else:
        logger.info("Reading particles from %s", inputParticlePath.resolve())
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection=inputParticles,
                orderedEvents=False,
            )
        )

    selectedParticles = "particles_selected"
    ptclSelector = ParticleSelector(
        level=acts.logging.INFO,
        inputParticles=inputParticles,
        outputParticles=selectedParticles,
        removeNeutral=True,
        absEtaMax=2.5,
        rhoMax=4.0 * u.mm,
        ptMin=500 * u.MeV,
    )
    s.addAlgorithm(ptclSelector)

    trackParameters = "trackparameters"
    if inputTrackSummary is None or inputParticlePath is None:
        logger.info("Using smeared particles")

        ptclSmearing = ParticleSmearing(
            level=acts.logging.INFO,
            inputParticles=selectedParticles,
            outputTrackParameters=trackParameters,
            randomNumbers=rnd,
        )
        s.addAlgorithm(ptclSmearing)
        associatedParticles = selectedParticles
    else:
        logger.info("Reading track summary from %s", inputTrackSummary.resolve())
        assert inputTrackSummary.exists()
        associatedParticles = "associatedTruthParticles"
        trackSummaryReader = RootTrajectorySummaryReader(
            level=acts.logging.VERBOSE,
            outputTracks="fittedTrackParameters",
            outputParticles=associatedParticles,
            filePath=str(inputTrackSummary.resolve()),
            orderedEvents=False,
        )
        s.addReader(trackSummaryReader)

        s.addAlgorithm(
            TrackSelector(
                level=acts.logging.INFO,
                inputTrackParameters=trackSummaryReader.config.outputTracks,
                outputTrackParameters=trackParameters,
                outputTrackIndices="outputTrackIndices",
                removeNeutral=True,
                absEtaMax=2.5,
                loc0Max=4.0 * u.mm,  # rho max
                ptMin=500 * u.MeV,
            )
        )

    logger.info("Using vertex finder: %s", vertexFinder.name)

    return addVertexFitting(
        s,
        field,
        outputDirRoot=outputDir if outputRoot else None,
        associatedParticles=associatedParticles,
        trackParameters=trackParameters,
        vertexFinder=vertexFinder,
    )


if "__main__" == __name__:
    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    inputTrackSummary = None
    for p in ("tracksummary_fitter.root", "tracksummary_ckf.root"):
        p = Path(p)
        if p.exists():
            inputTrackSummary = p
            break

    runVertexFitting(
        field,
        vertexFinder=VertexFinder.Truth,
        inputParticlePath=inputParticlePath,
        inputTrackSummary=inputTrackSummary,
        outputDir=Path.cwd(),
    ).run()
