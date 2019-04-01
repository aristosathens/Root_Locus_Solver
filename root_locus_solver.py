'''
    Aristos Athens
    Root Locus Solver

    Example:
        from root_locus_solver import root_locus
        poles = [ -4 + 2j, -4 - 2j, 1 ]
        zeros = [ -1 ]
        values = root_locus(poles = poles, zeros = zeros)
'''

import numpy as np                      # Use for finding roots, evaluationg polynomials, finding polynomial derivatives
from collections import Counter         # Use to count number of root occurrences

__all__ = [ "root_locus" ]              # Exported functions

# =============================================================================
# Constants

_ROUNDING_DECIMAL_PLACE = 3             # Allow for floating point rounding errors.

_POSITIVE = 1                           # Options for direction of change in K.
_NEGATIVE = 2

# =============================================================================
# Helper Functions

def _generate_polynomial_from_roots(roots):
    '''
        Given roots as numpy array, return polynomial as numpy array.
        polyfromroots() is an old function so it gives coefficients in
            reverse order (lowest degree first).
        Use np.flip() to get correct order (highest degree first).
    '''
    return np.flip(np.polynomial.polynomial.polyfromroots(roots))

def _wrap_angle(angle):
    '''
        Return angle value that is between [-360, 360].
        Handles lists and single values.
        Maintains angle sign.
        angle is in degrees.
    '''
    if isinstance(angle, (list, np.ndarray)):
        return np.array([_wrap_angle(a) for a in angle])
    else:
        if angle < 0:
            return (angle % 360) - 360
        else:
            return angle % 360

def _get_distance(point1, point2):
    '''
        Get the distance between two points (complex numbers).
        Can input any complex numbers (points, poles, zeros).
        Returns distance as float.
    '''
    delta_real = point1.real - point2.real
    delta_imag = point1.imag - point2.imag
    return np.sqrt(delta_real**2 + delta_imag**2)

def _get_angle(point1, point2):
    r'''
        Get the angle between two points (complex numbers).
            This is angle of point2 w.r.t. point1.
            |    p2    
            |      \
            |       \
            |        p1 ------ x
            |_______________
            This is the angle of line p2-p1 relative to line p1-x.
            In this drawing, angle will be positive.

        Can input any complex numbers (points, poles, zeros).
        Returns angle in degrees.
    '''
    delta_real = point2.real - point1.real
    delta_j = point2.imag - point1.imag
    angle = np.arctan2(delta_j, delta_real)
    return _wrap_angle(np.degrees(angle))

def _get_angles_sum(point, poles_or_zeros):
    '''
        Calculate the sum the angles between point and poles (or zeros).
    '''
    angles = [_get_angle(pole, point) for pole in poles_or_zeros]
    return np.sum(angles)

# =============================================================================
# Root Locus Functions 

def _points_on_root_locus(points, b_coefficients, a_coefficients, K_degree):
    '''
        Return the subset of points that is correctly on the locus.
        points is the set of values to check.

        If K_degree is positive, equation is:
            angle(K*b(s)/a(s)) = 180
        
        If K_degree is negative, equation is:
            angle(K*b(s)/a(s)) = 0
    '''
    # Recalculate poles and zeros from polynomials
    zeros = np.roots(b_coefficients)
    poles = np.roots(a_coefficients)

    sum_zero_angles = np.array([_get_angles_sum(point, zeros) for point in points])
    sum_pole_angles = np.array([_get_angles_sum(point, poles) for point in points])
    difference = sum_zero_angles - sum_pole_angles

    if K_degree == _POSITIVE:
        values = _wrap_angle(difference - 180)
    elif K_degree == _NEGATIVE:
        values = _wrap_angle(difference - 0)

    # Check for values are 0, within some rounding error
    valid_points = [v == 0 for v in np.round(values, decimals=_ROUNDING_DECIMAL_PLACE)]
    return points[valid_points]

def _angle_of_asymptote(n, m, l, K_degree):
    '''
        Find phi, the angle of the line of the asymptote.
        If K_degree is positive:
            phi = (180 - (l-1)360) / (n - m)
        If K_degree is negative:
            phi = 360(l - 1) / (n - m)
        Returns angle in degrees.
    '''
    if K_degree == _POSITIVE:
        return (180 + (l - 1)*360) / (n - m)
    elif K_degree == _NEGATIVE:
        return (l - 1)*360 / (n - m)

def _angles_of_all_asymptotes(poles, zeros, K_degree):
    '''
        Similar to _angle_of_asymptote, but compute all of the angles.
        Returns a numpy array of angles.
    '''
    n = len(poles)
    m = len(zeros)
    num_asymptotes = n - m

    angles = []
    for l in range(1, num_asymptotes+1):
        angles.append(_angle_of_asymptote(n, m, l, K_degree))
    
    return np.array(angles)

def _center_of_asymptotes(poles, zeros):
    '''
        Find alpha, the location from which asymptotes radiate. 
        Returns alpha as complex number.
    '''
    n = len(poles)
    m = len(zeros)
    return (np.sum(poles) - np.sum(zeros)) / (n - m)

def _angle_of_departure(poles, zeros, index, K_degree):
    '''
        Get angle of departure.
        If K_degree is positive:
            phi = 180 - (sum(t) - sum(phi))
        If K_degree is negative:
            phi = 0 - (sum(t) - sum(phi))
        Returns angle in degrees.
    '''
    # This is the point of interest
    p = poles[index]
    root_angles = [_get_angle(pole, p) for pole in poles if pole != p]
    zero_angles = [_get_angle(zero, p) for zero in zeros]
    if K_degree == _POSITIVE:
        return _wrap_angle(180 - (np.sum(zero_angles) - np.sum(root_angles)))
    elif K_degree == _NEGATIVE:
        return _wrap_angle(0 - (np.sum(zero_angles) - np.sum(root_angles)))

def _get_all_angles_of_departure(poles, zeros, K_degree):
    '''
        Similar to _angle_of_departure(), but does it for all poles.
        Returns numpy array of angles.
    '''
    angles = []
    for index in range(len(poles)):
        angles.append(_angle_of_departure(poles, zeros, index, K_degree))
    return angles

def _get_real_axis_points(b_coefficients, a_coefficients, K_degree):
    '''
        Rule 6.
        Get points where branches join/depart the real axis.
        Solve:
            b*a' - a*b' = 0
        Then check if points on root locus using angles
    '''
    d_b = np.polyder(b_coefficients)
    d_a = np.polyder(a_coefficients)
    eq = np.polysub(np.polymul(b_coefficients, d_a), np.polymul(a_coefficients, d_b))
    roots = np.roots(eq)

    # Remove any points not actually on root locus
    return np.real(_points_on_root_locus(roots, b_coefficients, a_coefficients, K_degree))
    
def _get_real_axis_angles(real_axis_points, b_coefficients, a_coefficients, K_degree):
    '''
        Rule 6.
        Find the angles of departure and approach for points on the real axis.
        Use equation:
            (180 + 360(l - 1)) / q, where q is number of roots for each real_axis_point
    '''
    angles = []
    for point in real_axis_points:
        q = _check_multiplicity(point, b_coefficients, a_coefficients, K_degree)
        angles.append([(180 + 360*(l - 1)) / q for l in range(q)])
    return angles

def _check_multiplicity(point, b_coefficients, a_coefficients, K_degree):
    '''
        Check how many branches approach/leave this point.
        If K_degree positive:
            K = -a(s)/b(s)
        If K_degree negative:
            K = a(s)/b(s)
    '''
    b_eval = np.polyval(b_coefficients, point)
    a_eval = np.polyval(a_coefficients, point)

    if K_degree == _POSITIVE:
        K = -a_eval / b_eval
    elif K_degree == _NEGATIVE:
        K = a_eval / b_eval

    # Find new roots of a(s) + Kb(s), for this value of K
    equation = np.polyadd(a_coefficients, np.polymul(K, b_coefficients))
    roots = np.roots(equation)

    # Only look at real value
    roots = [round(root.real, _ROUNDING_DECIMAL_PLACE) for root in roots]

    # Find max multipliciy of roots
    c = Counter(roots)
    c = {pair[0] : pair[1] for pair in c.items() if pair[1] > 1}
    return max(list(c.values()))

# =============================================================================
# Public Functions

def root_locus(b_coefficients = None, a_coefficients = None, zeros = None, poles = None, K_degree = "positive"):
    '''
        Get root locus info, given either coefficients for b(s), a(s) or the zeros/poles themselves.
        Must pass in b_coefficients and a_coefficients OR zeros and poles. No other combinations are allowed.
        K_degree is 'positive' or 'negative'. It determines direction of change in K.
        Returns dict of values.
    '''
    # If given cofficients, calculate zeros and poles.
    if b_coefficients == None and a_coefficients == None \
            and poles != None and zeros != None:
        b_coefficients = _generate_polynomial_from_roots(zeros)
        a_coefficients = _generate_polynomial_from_roots(poles)
    # If given poles and zeros, calculate coefficients.
    elif poles == None and zeros == None \
            and b_coefficients != None and a_coefficients != None:
        zeros = np.roots(b_coefficients)
        poles = np.roots(a_coefficients)
    else:
        raise Exception("Error: Passed in wrong combination. Requires either poles + zeros or a_coefficients + b_coefficients, not both.")

    if K_degree == "positive":
        K_degree = _POSITIVE
    elif K_degree == "negative":
        K_degree = _NEGATIVE
    else:
        raise Exception("Error: K_degree must be 'positive' or 'negative'.")

    asymptote_centroid = _center_of_asymptotes(poles, zeros)
    asymptote_angles = _angles_of_all_asymptotes(poles, zeros, K_degree)
    departure_angles = _get_all_angles_of_departure(poles, zeros, K_degree)
    real_axis_points = _get_real_axis_points(b_coefficients, a_coefficients, K_degree)
    real_axis_angles = _get_real_axis_angles(real_axis_points, b_coefficients, a_coefficients, K_degree)

    return {
        "poles" : poles,
        "zeros" : zeros,
        "asymptote_centroid" : asymptote_centroid,
        "asymptote_angles" : asymptote_angles,
        "departure_angles" : departure_angles,
        "real_axis_points" : real_axis_points,
        "real_axis_angles" : real_axis_angles,
    }
