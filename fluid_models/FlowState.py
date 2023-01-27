"""
Function:
Author: Luke Bartholomew
Edits:
"""
class FlowState():
    def __init__(self, state, vel_x = 0.0):
        self.fluid_state = state
        self.vel_x = vel_x