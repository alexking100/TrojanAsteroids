#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 21:13:05 2022

@author: alexking
"""
import numpy as np

myarray = np.array([[0,1,2,3],[4,5,6,7]])

np.savetxt('test.txt', myarray)