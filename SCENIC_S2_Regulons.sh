#!/bin/bash

export PATH=/BIGDATA1/pku_hkdeng_1/R/R-3.4.2/bin:$PATH

R CMD BATCH --no-save --no-restore S2_Regulons.R
