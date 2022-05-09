#!/bin/bash

set -e

outdir=$1
[ -z "$outdir" ] && outdir=physmon
mkdir -p $outdir

refdir=CI/physmon/reference
refcommit=$(cat $refdir/commit)
commit=$(git rev-parse --short HEAD)

echo "::group::Generate validation dataset"
CI/physmon/physmon.py $outdir 2>&1 > $outdir/run.log
echo "::endgroup::"

set +e

ec=0

function run() {
    a=$1
    b=$2

    echo "::group::Comparing $a vs. $b"

    histcmp \
        --label-reference=$refcommit \
        --label-monitored=$commit \
        "$@"

    ec=$(($ec | $?))

    echo "::endgroup::"
}


run \
    $outdir/performance_ckf_tracks_truth_smeared.root \
    $refdir/performance_ckf_tracks_truth_smeared.root \
    --title "CKF truth smeared" \
    -c CI/physmon/ckf_truth_smeared.yml \
    -o $outdir/ckf_truth_smeared.html \
    -p $outdir/ckf_truth_smeared_plots \

run \
    $outdir/performance_ckf_tracks_truth_estimated.root \
    $refdir/performance_ckf_tracks_truth_estimated.root \
    --title "CKF truth estimated" \
    -o $outdir/ckf_truth_estimated.html \
    -p $outdir/ckf_truth_estimated_plots \

run \
    $outdir/performance_ckf_tracks_seeded.root \
    $refdir/performance_ckf_tracks_seeded.root \
    --title "CKF seeded" \
    -o $outdir/ckf_seeded.html \
    -p $outdir/ckf_seeded_plots \

run \
    $outdir/performance_truth_tracking.root \
    $refdir/performance_truth_tracking.root \
    --title "Truth tracking" \
    -c CI/physmon/truth_tracking.yml \
    -o $outdir/truth_tracking.html \
    -p $outdir/truth_tracking_plots \


echo "::group::full_chain_odd"

thirdparty/OpenDataDetector/ci/full_chain_odd.py -o $outdir

build/bin/ActsAnalysisResidualsAndPulls \
    --predicted --filtered --smoothed --silent \
    -i $outdir/trackstates_ckf.root \
    -o $outdir/acts_analysis_residuals_and_pulls.root

histcmp \
    --title "full_chain_odd" \
    --label-reference=$refcommit \
    --label-monitored=$commit \
    $outdir/acts_analysis_residuals_and_pulls.root \
    CI/physmon/reference/acts_analysis_residuals_and_pulls.root

ec=$(($ec | $?))

echo "::endgroup::"


exit $ec
