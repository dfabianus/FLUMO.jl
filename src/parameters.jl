p_LDH = Dict(
    :a_n => 1.3343 ± 0.3837,
    :a_a => 12.0465 ± 1.3641,
    :b_n => -8.6824 ± 0.7182,
    :b_a => -16.7869 ± 2.5716,
    :p_u1 => [1.116222285164782e6, 0.6030700601335142],
    :p_u2 => [0.7903358250921652, 0.2, 0.12266210987282997],
    :beta1 => -0.25±0.03
)

p_GalOx = Dict(
    :a_n => 1.3343 ± 0.3837,
    :a_a => 12.0465 ± 0.3641,
    :b_n => -8.6824 ± 0.7182,
    :b_a => -16.7869 ± 0.5716,
    :a_ic => 0.0,
    :a_nc => 0.0,
    :a_cn => 1.0,
    :p_u1 => [1.116222285164782e6, 0.6030700601335142],
    :p_u2 => [0.7903358250921652, 0.2, 0.12266210987282997],
    :beta1 => -0.25±0.03
)

relative_errors(p,x) = p[1] * exp(-p[2]*x) + p[3]
relative_errors_native(x) = relative_errors([1.06, 17.66, 4.29],x)
#relative_errors_native(x) = relative_errors([32.3, 16.0, 1.54],x)
#relative_errors_aggregates(x) = relative_errors([35.38, 4.58, 9.43],x)
relative_errors_aggregates(x) = relative_errors([51.7, 29.0, 3.34],x)