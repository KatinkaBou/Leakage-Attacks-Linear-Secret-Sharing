# Simple Leakage Attacks on Secret Sharing Schemes

This GitHub repository contains the SageMath code in order to perform attacks against the leakage resilience of Shamir's and additive secret sharing schemes as presented in [BM26].
To run the code, simple execute the "attack.sage" code. It will run some pre-defined parameter sets. You can also run the attack with parameter sets of your choice.
Please make sure that your Sage version is up to date. We tested it on version 9.5 (Release Date: 2022-01-30) using Python 3.10.12.

#-- Example in Sage --#

    ..: load("attack.sage")

[BM26] Katharina Boudgoust, Mark Simkin. Scale, Round, Break: Simple Leakage Attacks on Secret Sharing Schemes
