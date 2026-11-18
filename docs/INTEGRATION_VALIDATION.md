# Official CleanLine integration validation (2026-11-18)

## Prior optimized-source validation

The source-package [validation report](https://github.com/MakotoMiyakoshi/CleanLine-GPT6-Optimized/blob/main/docs/VALIDATION_REPORT_20261118.md)
records exact equality (maximum absolute error 0, relative L2 error 0) of
cleaned samples and removed noise over 33 channels at 5–30 minutes. That
validation used the standalone plugin and a separate public CleanLine 2.00
reference. The source-package timing results remain in that report. They
were not remeasured for this official integration.

## New official-branch results

Ran `tests/test_cleanline_integration.m` with MATLAB R2024b Update 3 and
EEGLAB 2024.0 on deterministic 30-second, four-channel, 256-Hz synthetic
data. The baseline was this branch's unchanged official CleanLine 2.1 V2
(`pop_cleanline(...,'newversion',1)`). Results:

| Check | Result |
| --- | --- |
| Official V2 vs optimized Serial, single input | Exact; max absolute sample error 0 |
| Official V2 vs optimized Parallel, single input | Exact; max absolute sample error 0 |
| Official V2 vs optimized Serial, double input with channels 1 and 3 selected | Exact; max absolute sample error 0, using two valid official single-channel reference calls |
| Same nonconsecutive selection, optimized Parallel | Exact; max absolute sample error 0 |
| Original V2 integration route vs `pop_cleanline` | Exact |
| Historical V1 integration route vs `pop_cleanline` | Exact |
| Historical V1 ICA component route vs `pop_cleanline` | Exact |
| Generated EEGLAB command replay | Exact cleaned data |
| EEGLAB Tools registration | One CleanLine entry |

The test also confirmed that an unselected channel is unchanged. The
official V2 `parfor ch = chans` rejects nonconsecutive selections in MATLAB;
the optimized Parallel route handles them using consecutive loop indices.
This existing official V2 limitation is unchanged.

The adapted `tests/test_cleanline_gpt6_full_parity.m` is ready to compare
5- and 30-minute, 33-channel outputs directly against this checkout's
official V2. It was not run in this integration session because the
separately supplied validation `.set` dataset was unavailable here.
Interactive GUI layout and all parameter combinations remain for the
maintainer's alpha test.
