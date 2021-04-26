import ast
import numpy as np


def get_x_pos(pos):
    return parse_pos(pos)[0]


def get_y_pos(pos):
    return parse_pos(pos)[1]


def parse_pos(pos):
    try:
        return ast.literal_eval(pos)
    except (ValueError, SyntaxError):
        return [np.nan, np.nan]
