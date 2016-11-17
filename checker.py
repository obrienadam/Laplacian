#!/usr/bin/env python2

with open('result.txt', 'r') as f:
  result = f.read().replace('\n', '').split(',')
  result = [float(val) for val in result if val != ''] 
  print 'Max phi:', max(result)
  print 'Min phi:', min(result) 
