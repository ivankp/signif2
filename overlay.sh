#!/bin/bash

rxplot signif.root signif2.root -o signif.pdf \
  --yrange="0:5" --widths=4 -m '0.1:0.05:0.13:0.05' \
  --ticks-left --val-fmt="3.1f" --marker-size=2.2 --marker-color=1 \
  --xlabel-size=1.4 --xtitle-size=1.5 --xtitle-offset=1.0 \
  --ylabel-size=1.4 --ytitle-size=1.4 --ytitle-offset=1.0 \
  -r 'y/^.*/s\/#sqrt{s+b}' 'ng' \
  'nx/^(pT|m)_.*/& [GeV]' \
  'x/pT_([^ ]*)/p_{T}^{\1}' \
  'x/yAbs_yy/|y_{#gamma#gamma}|' \
  'x/cosTS_yy/|cos #theta_{#gamma#gamma}*|' \
  'x/N_j_(.*)/N_{ jets}^{ #geq\1 GeV}' \
  'x/Dphi_j_j/|#Delta#phi_{jj}|' \
  'x/m_jj/m_{jj}' \
  'x/Dy_j_j/|#Deltay_{jj}|'
