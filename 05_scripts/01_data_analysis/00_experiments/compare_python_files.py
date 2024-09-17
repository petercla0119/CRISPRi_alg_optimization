#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 09:57:50 2024

@author: Claire
"""

import difflib

with open("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_screenfiles/dual_guide_parser_v10-2.py") as file_1:
    file_1_parser = file_1.readlines()
    
with open("/Users/Claire/Library/CloudStorage/Box-Box/MSTP Research Rotations/Summer_2024_Might_Ward_Lab/CRISPRi_screenfiles/19_20_var_bp_dual_guide_parser_tool.py") as file_2:
        file_2_parser = file_2.readlines()
        
        
# Find and print the diff: 
diff = difflib.unified_diff(
    file_1_parser, file_2_parser, fromfile="dual_guide_parser_v10-2.py",
    tofile="19_20_var_bp_dual_guide_parser_tool.py", lineterm='')
# for line in difflib.unified_diff( 
#         file_1_parser, file_2_parser, fromfile='dual_guide_parser_v10-2.py',  
#         tofile='19_20_var_bp_dual_guide_parser_tool.py', lineterm=''): 
#     print(line) 
    
for lines in diff:
    print(lines)