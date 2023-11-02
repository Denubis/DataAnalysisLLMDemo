#!/usr/bin/env bash
# https://www.perplexity.ai/search/With-tar-I-7lTvFbXsQKiubqlGhSdSbQ?s=c

set -euo pipefail

rm -rf data
mkdir -p data

cat original_data/archive-name.tar.bz2.?? | tar -xjf - -C data

