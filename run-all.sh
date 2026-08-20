#!/bin/bash
Rscript quantile.R 2023-2024
Rscript data.R     2023-2024
DMCART_B=10000 Rscript dm-cart.R 2023-2024
Rscript visual.R   2023-2024
Rscript visual.R   2023-2024 lbw
