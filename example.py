#!/usr/bin/env python
from matplotlib.pyplot import show
import cartomap as cm

cm.plotCartoMap(states=False)
# states=False in case they haven't auto-downloaded detailed border data

show()
