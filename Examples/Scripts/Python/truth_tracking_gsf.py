#!/usr/bin/env python3
from pathlib import Path
from typing import Optional, Union

from acts.examples import Sequencer, GenericDetector, RootParticleReader

import acts

from acts import UnitConstants as u


def addGsfTracks(
    s: acts.examples.Sequencer,
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
):
    gsfOptions = {
        "maxComponents": 12,
        "abortOnError": False,
        "disableAllMaterialHandling": False,
    }

    gsfAlg = acts.examples.TrackFittingAlgorithm(
        level=acts.logging.INFO,
        inputMeasurements="measurements",
        inputSourceLinks="sourcelinks",
        inputProtoTracks="prototracks",
        inputInitialTrackParameters="estimatedparameters",
        outputTrajectories="gsf_trajectories",
        directNavigation=False,
        pickTrack=-1,
        trackingGeometry=trackingGeometry,
        fit=acts.examples.TrackFittingAlgorithm.makeGsfFitterFunction(
            trackingGeometry, field, **gsfOptions
        ),
    )

    s.addAlgorithm(gsfAlg)

    return s


def runGsfTracks(
    trackingGeometry,
    decorators,
    geometrySelection: Path,
    digiConfigFile: Path,
    field,
    outputDir: Path,
    truthSmearedSeeded=False,
    truthEstimatedSeeded=False,
    outputCsv=True,
    inputParticlePath: Optional[Path] = None,
    s=None,
):

    from particle_gun import addParticleGun, EtaConfig, PhiConfig, ParticleConfig
    from fatras import addFatras
    from digitization import addDigitization
    from seeding import addSeeding, SeedingAlgorithm

    s = s or acts.examples.Sequencer(
        events=100, numThreads=-1, logLevel=acts.logging.INFO
    )

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    if inputParticlePath is None:
        s = addParticleGun(
            s,
            EtaConfig(-2.0, 2.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, True),
            PhiConfig(0.0, 360.0 * u.degree),
            multiplicity=2,
            rnd=rnd,
        )
    else:
        acts.logging.getLogger("GSF Example").info(
            "Reading particles from %s", inputParticlePath.resolve()
        )
        assert inputParticlePath.exists()
        s.addReader(
            RootParticleReader(
                level=acts.logging.INFO,
                filePath=str(inputParticlePath.resolve()),
                particleCollection="particles_input",
                orderedEvents=False,
            )
        )

    s = addFatras(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
    )

    s = addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    s = addSeeding(
        s,
        trackingGeometry,
        field,
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
    )

    truthTrkFndAlg = acts.examples.TruthTrackFinder(
        level=acts.logging.INFO,
        inputParticles="truth_seeds_selected",
        inputMeasurementParticlesMap="measurement_particles_map",
        outputProtoTracks="prototracks",
    )

    s.addAlgorithm(truthTrkFndAlg)

    s = addGsfTracks(s, trackingGeometry, field)

    # Output
    s.addWriter(
        acts.examples.RootTrajectoryStatesWriter(
            level=acts.logging.INFO,
            inputTrajectories="gsf_trajectories",
            inputParticles="truth_seeds_selected",
            inputSimHits="simhits",
            inputMeasurementParticlesMap="measurement_particles_map",
            inputMeasurementSimHitsMap="measurement_simhits_map",
            filePath=str(outputDir / "trackstates_gsf.root"),
        )
    )

    s.addWriter(
        acts.examples.RootTrajectorySummaryWriter(
            level=acts.logging.INFO,
            inputTrajectories="gsf_trajectories",
            inputParticles="truth_seeds_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            filePath=str(outputDir / "tracksummary_gsf.root"),
        )
    )

    return s


if "__main__" == __name__:
    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    detector, trackingGeometry, decorators = GenericDetector.create()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    inputParticlePath = Path("particles.root")
    if not inputParticlePath.exists():
        inputParticlePath = None

    runGsfTracks(
        trackingGeometry,
        decorators,
        field=field,
        geometrySelection=srcdir
        / "Examples/Algorithms/TrackFinding/share/geoSelection-genericDetector.json",
        digiConfigFile=srcdir
        / "Examples/Algorithms/Digitization/share/default-smearing-config-generic.json",
        outputCsv=True,
        truthSmearedSeeded=False,
        truthEstimatedSeeded=False,
        inputParticlePath=inputParticlePath,
        outputDir=Path.cwd(),
    ).run()
