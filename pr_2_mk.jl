#INTERPOLACJA

## suma wart
function rozwiazanie_2(;
    m::Vector{Float64} = [-4.6, -4.5943, -4.5886, -4.5829, -4.5772, -4.5715, -4.5658, -4.5601, -4.5544, -4.5487, -4.543, -4.5373, -4.5316, -4.5259, -4.5202, -4.5145, -4.5088, -4.5031, -4.4974, -4.4917, -4.486, -4.4803, -4.4746, -4.4689, -4.4632, -4.4575, -4.4518, -4.4461, -4.4404, -4.4347, -4.429, -4.4233, -4.4176, -4.4119, -4.4062, -4.4005, -4.3948, -4.3891, -4.3834, -4.3777, -4.372, -4.3663, -4.3606, -4.3549, -4.3492, -4.3435, -4.3378, -4.3321, -4.3264, -4.3207, -4.315, -4.3093, -4.3036, -4.2979, -4.2922, -4.2865, -4.2808, -4.2751, -4.2694, -4.2637, -4.258, -4.2523, -4.2466, -4.2409, -4.2352, -4.2295, -4.2238, -4.2181, -4.2124, -4.2067, -4.201, -4.1953, -4.1896, -4.1839, -4.1782, -4.1725, -4.1668, -4.1611, -4.1554, -4.1497, -4.144, -4.1383, -4.1326, -4.1269, -4.1212, -4.1155, -4.1098, -4.1041, -4.0984, -4.0927],
    s::Vector{Float64} = [0.9181, 0.0832, 0.9099, 0.0934, 0.3633, 0.0321, 0.7866, 0.631, 0.5492, 0.3496, 0.2652, 0.6898, 0.264, 0.9838, 0.7483, 0.6867, 0.4623, 0.8455, 0.1527, 0.9919, 0.0016, 0.6075, 0.008, 0.1845, 0.3393, 0.4392, 0.5838, 0.4931, 0.0014, 0.7951, 0.092, 0.7579, 0.2956, 0.9166, 0.917, 0.3376, 0.0481, 0.9035, 0.6418, 0.6664, 0.603, 0.1425, 0.5517, 0.6841, 0.5058, 0.4153, 0.018, 0.255, 0.3559, 0.5227, 0.865, 0.0313, 0.73, 0.1262, 0.2263, 0.5154, 0.4437, 0.3205, 0.1881, 0.1671, 0.3574, 0.4303, 0.6071, 0.7636, 0.1639, 0.1147, 0.8315, 0.4347, 0.0426, 0.4521, 0.1507, 0.2625, 0.0196, 0.9987, 0.4847, 0.5893, 0.7616, 0.6314, 0.5737, 0.6568, 0.3523, 0.1118, 0.9396, 0.4268, 0.7448, 0.7097, 0.5352, 0.6985, 0.7222, 0.5939],
    t::Vector{Float64} = [-4.5772, -4.50082, -4.18105, -4.09783, -4.36801, -4.16737, -4.45921, -4.42558, -4.47289, -4.55041],
)
    #4.723087822917896
    out = zeros(length(t))
    T = m[2]-m[1]
    for i in 1:length(t)
        for n in 1:length(s)
            out[i] += sinc((t[i]-m[n])/T)*s[n]
        end
    end
    return sum(out)
end
rozwiazanie_2()

## suma wart
function rozwiazanie_2(;
    m::Vector{Float64} = [3.0, 3.0079, 3.0158, 3.0237, 3.0316, 3.0395, 3.0474, 3.0553, 3.0632, 3.0711, 3.079, 3.0869, 3.0948, 3.1027, 3.1106, 3.1185, 3.1264, 3.1343, 3.1422, 3.1501, 3.158, 3.1659, 3.1738, 3.1817, 3.1896, 3.1975, 3.2054, 3.2133, 3.2212, 3.2291, 3.237, 3.2449, 3.2528, 3.2607, 3.2686, 3.2765, 3.2844, 3.2923, 3.3002, 3.3081, 3.316, 3.3239, 3.3318, 3.3397, 3.3476, 3.3555, 3.3634, 3.3713, 3.3792, 3.3871, 3.395, 3.4029, 3.4108, 3.4187, 3.4266, 3.4345, 3.4424, 3.4503, 3.4582, 3.4661, 3.474, 3.4819, 3.4898, 3.4977, 3.5056, 3.5135, 3.5214, 3.5293, 3.5372, 3.5451, 3.553, 3.5609, 3.5688, 3.5767, 3.5846, 3.5925, 3.6004, 3.6083, 3.6162, 3.6241, 3.632, 3.6399, 3.6478, 3.6557, 3.6636, 3.6715, 3.6794, 3.6873, 3.6952, 3.7031, 3.711, 3.7189, 3.7268, 3.7347, 3.7426, 3.7505, 3.7584],
    s::Vector{Float64} = [0.4248, 0.1917, 0.435, 0.6436, 0.3908, 0.9747, 0.7978, 0.3622, 0.1255, 0.5443, 0.0586, 0.3332, 0.8076, 0.9643, 0.8925, 0.2312, 0.2874, 0.5902, 0.5251, 0.6284, 0.5546, 0.3058, 0.1847, 0.9289, 0.0629, 0.8591, 0.946, 0.1711, 0.8057, 0.3376, 0.2828, 0.0719, 0.2345, 0.2813, 0.6257, 0.6766, 0.1322, 0.2776, 0.8441, 0.3391, 0.2295, 0.3125, 0.1133, 0.234, 0.7742, 0.0954, 0.6206, 0.3865, 0.4015, 0.0791, 0.2481, 0.6652, 0.0528, 0.8686, 0.8018, 0.7351, 0.0872, 0.3354, 0.1501, 0.7375, 0.3727, 0.6676, 0.6272, 0.9821, 0.7076, 0.1242, 0.7929, 0.9128, 0.5809, 0.178, 0.2446, 0.3586, 0.9424, 0.3758, 0.7151, 0.0623, 0.9864, 0.1462, 0.7532, 0.0097, 0.1662, 0.13, 0.4269, 0.0796, 0.3998, 0.1178, 0.0847, 0.7395, 0.4215, 0.4792, 0.9564, 0.0066, 0.4103, 0.056, 0.7725, 0.1567, 0.6203],
    t::Vector{Float64} = [3.16906, 3.03792, 3.38236, 3.07821, 3.06478, 3.38947, 3.36735, 3.62252, 3.04029, 3.29625, 3.49296],
)
    #4.776192451258249
    out = zeros(length(t))
    T = m[2] - m[1]
    for i in 1:length(t)
        for n in 1:length(s)
            out[i] += s[n]*sinc((t[i]-m[n])/T)
        end
    end
    return sum(out)
end
rozwiazanie_2()

## suma wart
function rozwiazanie_2(;
    m::Vector{Float64} = [-0.6, -0.5929, -0.5858, -0.5787, -0.5716, -0.5645, -0.5574, -0.5503, -0.5432, -0.5361, -0.529, -0.5219, -0.5148, -0.5077, -0.5006, -0.4935, -0.4864, -0.4793, -0.4722, -0.4651, -0.458, -0.4509, -0.4438, -0.4367, -0.4296, -0.4225, -0.4154, -0.4083, -0.4012, -0.3941, -0.387, -0.3799, -0.3728, -0.3657, -0.3586, -0.3515, -0.3444, -0.3373, -0.3302, -0.3231, -0.316, -0.3089, -0.3018, -0.2947, -0.2876, -0.2805, -0.2734, -0.2663, -0.2592, -0.2521, -0.245, -0.2379, -0.2308, -0.2237, -0.2166, -0.2095, -0.2024, -0.1953, -0.1882, -0.1811, -0.174, -0.1669, -0.1598, -0.1527, -0.1456, -0.1385, -0.1314, -0.1243, -0.1172],
    s::Vector{Float64} = [0.2858, 0.0776, 0.1794, 0.8403, 0.0362, 0.4326, 0.8338, 0.5308, 0.0357, 0.5726, 0.0807, 0.4245, 0.4592, 0.8551, 0.5058, 0.6087, 0.6029, 0.9367, 0.9896, 0.6597, 0.0684, 0.6888, 0.3553, 0.3369, 0.9765, 0.4928, 0.0022, 0.8559, 0.6615, 0.4097, 0.5143, 0.2674, 0.0965, 0.5254, 0.3005, 0.4864, 0.5025, 0.5806, 0.4303, 0.9501, 0.7656, 0.259, 0.9919, 0.915, 0.0163, 0.3192, 0.094, 0.102, 0.2252, 0.1916, 0.4152, 0.8388, 0.9572, 0.2236, 0.7749, 0.14, 0.7385, 0.5356, 0.6622, 0.382, 0.6871, 0.8261, 0.9842, 0.703, 0.8362, 0.4976, 0.7063, 0.1424, 0.1214],
    t::Vector{Float64} = [-0.36144, -0.41185, -0.50557, -0.2237, -0.32665, -0.52261, -0.43386, -0.34298, -0.5858, -0.22015, -0.41824, -0.41966],
)
    #4.993366366477644
    out = zeros(length(t))
    T = m[2] - m[1]
    for i in 1:length(t)
        for n in 1:length(s)
            out[i] += s[n]*sinc((t[i]-m[n])/T)
        end
    end
    return sum(out)
end
rozwiazanie_2()

## suma wart
function rozwiazanie_2(;
    m::Vector{Float64} = [4.5, 4.5073, 4.5146, 4.5219, 4.5292, 4.5365, 4.5438, 4.5511, 4.5584, 4.5657, 4.573, 4.5803, 4.5876, 4.5949, 4.6022, 4.6095, 4.6168, 4.6241, 4.6314, 4.6387, 4.646, 4.6533, 4.6606, 4.6679, 4.6752, 4.6825, 4.6898, 4.6971, 4.7044, 4.7117, 4.719, 4.7263, 4.7336, 4.7409, 4.7482, 4.7555, 4.7628, 4.7701, 4.7774, 4.7847, 4.792, 4.7993, 4.8066, 4.8139, 4.8212, 4.8285, 4.8358, 4.8431, 4.8504, 4.8577, 4.865, 4.8723, 4.8796, 4.8869, 4.8942, 4.9015, 4.9088, 4.9161, 4.9234, 4.9307, 4.938, 4.9453, 4.9526, 4.9599, 4.9672, 4.9745, 4.9818, 4.9891, 4.9964, 5.0037, 5.011, 5.0183, 5.0256, 5.0329, 5.0402, 5.0475, 5.0548, 5.0621, 5.0694, 5.0767, 5.084, 5.0913, 5.0986, 5.1059, 5.1132, 5.1205, 5.1278, 5.1351, 5.1424],
    s::Vector{Float64} = [0.8663, 0.2915, 0.0418, 0.6486, 0.6286, 0.8537, 0.534, 0.6275, 0.3427, 0.5091, 0.2794, 0.8473, 0.9681, 0.52, 0.7633, 0.1209, 0.3497, 0.0504, 0.1535, 0.2793, 0.9164, 0.4056, 0.5923, 0.4826, 0.765, 0.056, 0.004, 0.7672, 0.3431, 0.7192, 0.1289, 0.1212, 0.8814, 0.7124, 0.3345, 0.5755, 0.2649, 0.2322, 0.4161, 0.8611, 0.3825, 0.2227, 0.7007, 0.9048, 0.6419, 0.6837, 0.7922, 0.7036, 0.7466, 0.844, 0.6692, 0.8124, 0.09, 0.7421, 0.7172, 0.2761, 0.9028, 0.4343, 0.8197, 0.3989, 0.5002, 0.971, 0.0004, 0.8821, 0.4301, 0.828, 0.7232, 0.4802, 0.2385, 0.9299, 0.0418, 0.2475, 0.0056, 0.9192, 0.7687, 0.0733, 0.3733, 0.7668, 0.0616, 0.7759, 0.7295, 0.1857, 0.3277, 0.0548, 0.4517, 0.251, 0.6162, 0.4389, 0.5913],
    t::Vector{Float64} = [4.69783, 4.94603, 4.64235, 4.72849, 4.94749, 4.82047, 4.61534, 4.87084, 4.72046, 5.04458, 4.74966, 4.89201, 4.71462],
)
    #7.021873602057676
    out = zeros(length(t))
    T = m[2]-m[1]
    for i in 1:length(t)
        for n in 1:length(s)
            out[i] += s[n]*sinc((t[i]-m[n])/T)
        end
    end
    return sum(out)
end
rozwiazanie_2()
