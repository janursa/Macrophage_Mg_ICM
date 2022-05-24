#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 24 18:11:03 2022

@author: matin
"""


import tellurium as te
import sys
import os
from pathlib import Path
dir_file = Path(__file__).resolve().parent
main_dir = os.path.join(dir_file,'..')
sys.path.insert(0,main_dir)
import json
sys.path.insert(0,main_dir)
from data.observations import observations,t2m,select_obs,packages
from models.params import fixed_params
from tools import dirs, tools
# from models.models import Macrophage
import matplotlib
import matplotlib.pyplot as plt
import tellurium as te
from models.models import Macrophage

# tools.activation_LPS()

obj = Macrophage('M1')
obj.run(params={},studies = select_obs(['eq_mg']))
print('completed')