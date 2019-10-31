#!/usr/bin/env python3

import GiniPackage as gp

a = gp.GiniObject()
a.compile("Data/sampledata.csv")
a.setTSNE2("Data/sample_tsne.csv")
a.setClusters("Data/Sample_Clusters.csv")
a.DrawGini()