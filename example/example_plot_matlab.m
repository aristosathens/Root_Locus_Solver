poles = [ -4 + 2j, -4 - 2j, 1 ];
zeros = [ -1 ];

b_coefficients = poly(zeros);
a_coefficients = poly(poles);

H = tf(b_coefficients,a_coefficients);
rlocus(H)