#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 21:13:05 2022

@author: alexking
"""
import numpy as np

arr = np.array([
	[0,0],
	[1,1],
	[2,2],
	[3,3],
	[4,4]
	])
norms = np.linalg.norm(arr, axis = 1)
print(arr)
print(norms)

print(np.amax(norms))