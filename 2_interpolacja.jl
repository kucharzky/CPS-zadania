## suma wartosci 

function rozwiazanie(;
    m::Vector{Float64} = [1.8, 1.8001, 1.8002, 1.8003, 1.8004, 1.8005, 1.8006, 1.8007, 1.8008, 1.8009, 1.801, 1.8011, 1.8012, 1.8013, 1.8014, 1.8015, 1.8016, 1.8017, 1.8018, 1.8019, 1.802, 1.8021, 1.8022, 1.8023, 1.8024, 1.8025, 1.8026, 1.8027, 1.8028, 1.8029, 1.803, 1.8031, 1.8032, 1.8033, 1.8034, 1.8035, 1.8036, 1.8037, 1.8038, 1.8039, 1.804, 1.8041, 1.8042, 1.8043, 1.8044, 1.8045, 1.8046, 1.8047, 1.8048, 1.8049, 1.805, 1.8051, 1.8052, 1.8053, 1.8054],
    s::Vector{Float64} = [0.0597, 0.8512, 0.8973, 0.0039, 0.8853, 0.4719, 0.622, 0.3112, 0.008, 0.1859, 0.7256, 0.1486, 0.6359, 0.5802, 0.567, 0.6844, 0.9525, 0.3516, 0.6431, 0.6968, 0.6276, 0.2729, 0.1153, 0.6208, 0.3319, 0.402, 0.4994, 0.6306, 0.746, 0.4849, 0.1749, 0.0994, 0.3751, 0.5838, 0.5855, 0.5086, 0.8393, 0.0931, 0.4736, 0.2923, 0.5472, 0.5224, 0.462, 0.935, 0.3491, 0.3401, 0.6853, 0.1801, 0.7504, 0.9484, 0.9799, 0.583, 0.6174, 0.4253, 0.1541],
    t::Vector{Float64} = [1.80232, 1.80035, 1.80221, 1.80379, 1.80534, 1.80425, 1.80113, 1.80214, 1.80226, 1.80131, 1.80449, 1.8023, 1.80304, 1.8036, 1.80026, 1.80437],
)
    # 6.48141032359112
    t_out = zeros(length(t))
    T = m[2]-m[1]
    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i]+=sinc((t[i]-m[n])/T)*s[n]
        end
    end
    return sum(t_out)
end

rozwiazanie()

## suma wartosci
          
function rozwiazanie(;
    m::Vector{Float64} = [4.4, 4.4045, 4.409, 4.4135, 4.418, 4.4225, 4.427, 4.4315, 4.436, 4.4405, 4.445, 4.4495, 4.454, 4.4585, 4.463, 4.4675, 4.472, 4.4765, 4.481, 4.4855, 4.49, 4.4945, 4.499, 4.5035, 4.508, 4.5125, 4.517, 4.5215, 4.526, 4.5305, 4.535, 4.5395, 4.544, 4.5485, 4.553, 4.5575, 4.562, 4.5665, 4.571, 4.5755, 4.58, 4.5845, 4.589, 4.5935, 4.598, 4.6025, 4.607, 4.6115, 4.616, 4.6205, 4.625, 4.6295, 4.634, 4.6385, 4.643, 4.6475, 4.652, 4.6565, 4.661, 4.6655, 4.67, 4.6745, 4.679, 4.6835, 4.688, 4.6925, 4.697, 4.7015, 4.706, 4.7105, 4.715, 4.7195],
    s::Vector{Float64} = [0.7841, 0.8248, 0.58, 0.6996, 0.0065, 0.6845, 0.9185, 0.1753, 0.2018, 0.6241, 0.3312, 0.7459, 0.5527, 0.8732, 0.6315, 0.3305, 0.1869, 0.9963, 0.4861, 0.7649, 0.5217, 0.8049, 0.1221, 0.5738, 0.7062, 0.5835, 0.3065, 0.1142, 0.657, 0.1897, 0.6586, 0.764, 0.6832, 0.1583, 0.2515, 0.5384, 0.4849, 0.2, 0.8091, 0.5257, 0.3149, 0.2359, 0.358, 0.3318, 0.2422, 0.2952, 0.2466, 0.1261, 0.1148, 0.3602, 0.512, 0.4256, 0.9561, 0.7895, 0.1664, 0.2788, 0.3812, 0.4581, 0.9192, 0.8288, 0.8774, 0.3023, 0.2171, 0.6666, 0.2019, 0.0369, 0.3516, 0.6421, 0.698, 0.1454, 0.46, 0.3891],
    t::Vector{Float64} = [4.51475, 4.4819, 4.52645, 4.6925, 4.63715, 4.44095, 4.60655, 4.6142, 4.46255],
)
    # 4.336861844050346
    t_out = zeros(length(t))
    T = m[2] - m[1]
    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i] += sinc((t[i]-m[n])/T)*s[n]
        end
    end

    return sum(t_out)
end

rozwiazanie()

## suma wartosci
          
function rozwiazanie(;
    m::Vector{Float64} = [4.1, 4.1072, 4.1144, 4.1216, 4.1288, 4.136, 4.1432, 4.1504, 4.1576, 4.1648, 4.172, 4.1792, 4.1864, 4.1936, 4.2008, 4.208, 4.2152, 4.2224, 4.2296, 4.2368, 4.244, 4.2512, 4.2584, 4.2656, 4.2728, 4.28, 4.2872, 4.2944, 4.3016, 4.3088, 4.316, 4.3232, 4.3304, 4.3376, 4.3448, 4.352, 4.3592, 4.3664, 4.3736, 4.3808, 4.388, 4.3952, 4.4024, 4.4096, 4.4168, 4.424, 4.4312, 4.4384, 4.4456, 4.4528, 4.46, 4.4672, 4.4744, 4.4816, 4.4888, 4.496, 4.5032, 4.5104, 4.5176, 4.5248, 4.532, 4.5392, 4.5464, 4.5536, 4.5608, 4.568, 4.5752, 4.5824, 4.5896, 4.5968, 4.604, 4.6112, 4.6184, 4.6256, 4.6328, 4.64, 4.6472, 4.6544, 4.6616, 4.6688, 4.676, 4.6832, 4.6904, 4.6976, 4.7048],
    s::Vector{Float64} = [0.376, 0.2085, 0.3809, 0.9032, 0.1743, 0.4908, 0.3514, 0.8997, 0.8108, 0.5514, 0.6091, 0.134, 0.0724, 0.1042, 0.4998, 0.6104, 0.8204, 0.0194, 0.2213, 0.8509, 0.454, 0.8178, 0.1795, 0.0518, 0.1221, 0.4801, 0.9229, 0.3352, 0.8223, 0.9667, 0.7864, 0.6938, 0.6803, 0.846, 0.6683, 0.5426, 0.5906, 0.0754, 0.5072, 0.3286, 0.7547, 0.1261, 0.4, 0.2293, 0.0172, 0.6388, 0.0012, 0.9548, 0.5979, 0.995, 0.9839, 0.9722, 0.3481, 0.4941, 0.8777, 0.2267, 0.5505, 0.739, 0.2058, 0.6466, 0.6428, 0.1309, 0.1771, 0.0578, 0.6268, 0.8046, 0.731, 0.4474, 0.7406, 0.1344, 0.2377, 0.3806, 0.8368, 0.3124, 0.1126, 0.3825, 0.7853, 0.73, 0.1042, 0.3225, 0.1509, 0.0499, 0.1424, 0.8731, 0.1251],
    t::Vector{Float64} = [4.22384, 4.52192, 4.44344],
)
    # 0.9937904266412192
    T = m[2] - m[1]
    t_out = zeros(Float64,length(t))

    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i]+=sinc((t[i]-m[n])/T)*s[n]
        end
    end
    return sum(t_out)
end

rozwiazanie()

## suma wartosci
          
function rozwiazanie(;
    m::Vector{Float64} = [3.1, 3.1021, 3.1042, 3.1063, 3.1084, 3.1105, 3.1126, 3.1147, 3.1168, 3.1189, 3.121, 3.1231, 3.1252, 3.1273, 3.1294, 3.1315, 3.1336, 3.1357, 3.1378, 3.1399, 3.142, 3.1441, 3.1462, 3.1483, 3.1504, 3.1525, 3.1546, 3.1567, 3.1588, 3.1609, 3.163, 3.1651, 3.1672, 3.1693, 3.1714, 3.1735, 3.1756, 3.1777, 3.1798, 3.1819, 3.184, 3.1861, 3.1882, 3.1903, 3.1924, 3.1945, 3.1966, 3.1987, 3.2008, 3.2029, 3.205, 3.2071, 3.2092, 3.2113, 3.2134, 3.2155],
    s::Vector{Float64} = [0.8032, 0.4027, 0.1633, 0.016, 0.3568, 0.5183, 0.8469, 0.6265, 0.0171, 0.2286, 0.5351, 0.2018, 0.3287, 0.3647, 0.6046, 0.885, 0.1782, 0.6175, 0.1859, 0.2014, 0.0054, 0.4564, 0.8453, 0.5645, 0.5313, 0.8923, 0.4677, 0.5117, 0.7014, 0.7307, 0.1907, 0.3084, 0.2579, 0.8727, 0.9254, 0.2269, 0.4041, 0.3393, 0.2581, 0.9764, 0.5018, 0.3339, 0.0475, 0.9538, 0.8738, 0.1682, 0.2405, 0.0006, 0.1616, 0.2661, 0.0293, 0.3019, 0.1409, 0.6696, 0.2294, 0.6789],
    t::Vector{Float64} = [3.1294, 3.17539, 3.12079, 3.14032, 3.13171, 3.16174],
)
    # 2.9447803833684785
    T = m[2] - m[1]
    t_out = zeros(Float64,length(t))
    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i]+=sinc((t[i]-m[n])/T)*s[n]
        end
    end
    return sum(t_out)
end

rozwiazanie()

## suma wartosci

          
function rozwiazanie(;
    m::Vector{Float64} = [-4.8, -4.7959, -4.7918, -4.7877, -4.7836, -4.7795, -4.7754, -4.7713, -4.7672, -4.7631, -4.759, -4.7549, -4.7508, -4.7467, -4.7426, -4.7385, -4.7344, -4.7303, -4.7262, -4.7221, -4.718, -4.7139, -4.7098, -4.7057, -4.7016, -4.6975, -4.6934, -4.6893, -4.6852, -4.6811, -4.677, -4.6729, -4.6688, -4.6647, -4.6606, -4.6565, -4.6524, -4.6483, -4.6442, -4.6401, -4.636, -4.6319, -4.6278, -4.6237, -4.6196, -4.6155, -4.6114, -4.6073, -4.6032, -4.5991, -4.595, -4.5909, -4.5868, -4.5827, -4.5786, -4.5745, -4.5704, -4.5663, -4.5622, -4.5581, -4.554, -4.5499, -4.5458, -4.5417, -4.5376, -4.5335],
    s::Vector{Float64} = [0.5913, 0.0411, 0.8095, 0.2103, 0.922, 0.1429, 0.368, 0.5462, 0.1496, 0.0985, 0.1344, 0.9982, 0.573, 0.8254, 0.1913, 0.949, 0.7859, 0.8826, 0.2834, 0.8739, 0.9213, 0.0294, 0.8099, 0.116, 0.5774, 0.3571, 0.1756, 0.4913, 0.4748, 0.5672, 0.2781, 0.8256, 0.619, 0.2865, 0.4772, 0.4735, 0.3861, 0.6674, 0.5678, 0.6751, 0.8, 0.4981, 0.5631, 0.7581, 0.7379, 0.6619, 0.4964, 0.775, 0.0875, 0.9499, 0.9315, 0.7798, 0.5945, 0.4218, 0.3751, 0.8099, 0.1278, 0.1907, 0.2561, 0.7943, 0.7429, 0.962, 0.576, 0.2439, 0.6458, 0.69],
    t::Vector{Float64} = [-4.77417, -4.56876, -4.78237, -4.60935, -4.57901, -4.77376],
)
    # 3.622401658662879
    T = m[2] - m[1]
    t_out = zeros(Float64,length(t))

    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i]+= sinc((t[i]-m[n])/T)*s[n]
        end
    end

    return sum(t_out)
end

rozwiazanie()

## suma wartosci

          
function rozwiazanie(;
    m::Vector{Float64} = [-2.2, -2.1936, -2.1872, -2.1808, -2.1744, -2.168, -2.1616, -2.1552, -2.1488, -2.1424, -2.136, -2.1296, -2.1232, -2.1168, -2.1104, -2.104, -2.0976, -2.0912, -2.0848, -2.0784, -2.072, -2.0656, -2.0592, -2.0528, -2.0464, -2.04, -2.0336, -2.0272, -2.0208, -2.0144, -2.008, -2.0016, -1.9952, -1.9888, -1.9824, -1.976, -1.9696, -1.9632, -1.9568, -1.9504, -1.944, -1.9376, -1.9312, -1.9248, -1.9184, -1.912, -1.9056, -1.8992, -1.8928, -1.8864, -1.88, -1.8736, -1.8672, -1.8608, -1.8544, -1.848, -1.8416, -1.8352, -1.8288, -1.8224, -1.816, -1.8096, -1.8032, -1.7968, -1.7904, -1.784, -1.7776, -1.7712, -1.7648, -1.7584, -1.752, -1.7456, -1.7392, -1.7328, -1.7264, -1.72, -1.7136, -1.7072, -1.7008, -1.6944, -1.688, -1.6816, -1.6752, -1.6688, -1.6624, -1.656, -1.6496, -1.6432, -1.6368],
    s::Vector{Float64} = [0.4608, 0.4796, 0.7411, 0.9659, 0.7751, 0.1624, 0.3913, 0.4888, 0.9113, 0.0129, 0.9442, 0.698, 0.7955, 0.6285, 0.8501, 0.8154, 0.2151, 0.7949, 0.0427, 0.4573, 0.9127, 0.0821, 0.6835, 0.3361, 0.6968, 0.6493, 0.1388, 0.2756, 0.0251, 0.8027, 0.9432, 0.5072, 0.1192, 0.3725, 0.7293, 0.58, 0.1931, 0.7929, 0.9802, 0.9441, 0.7254, 0.5431, 0.2492, 0.4598, 0.705, 0.3196, 0.3963, 0.0943, 0.7817, 0.0271, 0.8379, 0.6201, 0.265, 0.1583, 0.7979, 0.8926, 0.6858, 0.673, 0.0289, 0.3056, 0.4629, 0.9885, 0.7715, 0.9151, 0.0263, 0.2515, 0.8432, 0.7596, 0.5007, 0.8033, 0.6818, 0.721, 0.3005, 0.2888, 0.0908, 0.5942, 0.7545, 0.7923, 0.9862, 0.691, 0.0039, 0.838, 0.2256, 0.233, 0.6945, 0.4373, 0.9173, 0.3584, 0.1759],
    t::Vector{Float64} = [-2.00288, -2.19168, -1.69184, -1.82688, -1.84224, -1.64832, -1.84352, -1.7072, -1.73984],
)
    # 4.799202605542257
    t_out = zeros(Float64,length(t))
end

## suma wartosci

          
function rozwiazanie(;
    m::Vector{Float64} = [-0.9, -0.8973, -0.8946, -0.8919, -0.8892, -0.8865, -0.8838, -0.8811, -0.8784, -0.8757, -0.873, -0.8703, -0.8676, -0.8649, -0.8622, -0.8595, -0.8568, -0.8541, -0.8514, -0.8487, -0.846, -0.8433, -0.8406, -0.8379, -0.8352, -0.8325, -0.8298, -0.8271, -0.8244, -0.8217, -0.819, -0.8163, -0.8136, -0.8109, -0.8082, -0.8055, -0.8028, -0.8001, -0.7974, -0.7947, -0.792, -0.7893, -0.7866, -0.7839, -0.7812, -0.7785, -0.7758, -0.7731, -0.7704, -0.7677, -0.765, -0.7623, -0.7596, -0.7569, -0.7542, -0.7515, -0.7488, -0.7461, -0.7434, -0.7407, -0.738, -0.7353, -0.7326, -0.7299, -0.7272, -0.7245, -0.7218, -0.7191, -0.7164, -0.7137, -0.711, -0.7083, -0.7056, -0.7029, -0.7002, -0.6975, -0.6948, -0.6921, -0.6894, -0.6867, -0.684],
    s::Vector{Float64} = [0.6403, 0.6013, 0.8024, 0.9689, 0.8877, 0.7264, 0.1324, 0.9595, 0.2465, 0.4998, 0.1688, 0.1266, 0.4028, 0.8224, 0.689, 0.2521, 0.8128, 0.9593, 0.3154, 0.2478, 0.5553, 0.7324, 0.2829, 0.5358, 0.7789, 0.9217, 0.7037, 0.7497, 0.3548, 0.4591, 0.2717, 0.1374, 0.5683, 0.8673, 0.7903, 0.0367, 0.577, 0.7836, 0.6494, 0.2879, 0.2323, 0.1425, 0.1869, 0.9871, 0.5738, 0.3331, 0.8562, 0.8163, 0.5119, 0.3581, 0.6524, 0.0072, 0.1432, 0.1106, 0.7423, 0.3598, 0.1451, 0.0594, 0.6144, 0.6956, 0.8862, 0.6086, 0.526, 0.9395, 0.2394, 0.2636, 0.4926, 0.0617, 0.2644, 0.1865, 0.6351, 0.3529, 0.9189, 0.1851, 0.1673, 0.0389, 0.565, 0.6436, 0.3543, 0.7413, 0.1335],
    t::Vector{Float64} = [-0.70398, -0.74772, -0.77742],
)
    # 1.3030357632665153
    T = m[2]-m[1]
    t_out = zeros(length(t))

    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i]+=sinc((t[i]-m[n])/T)*s[n]
        end
    end

    return sum(t_out) 
end

rozwiazanie()

## suma wartosci

          
function rozwiazanie(;
    m::Vector{Float64} = [2.7, 2.7011, 2.7022, 2.7033, 2.7044, 2.7055, 2.7066, 2.7077, 2.7088, 2.7099, 2.711, 2.7121, 2.7132, 2.7143, 2.7154, 2.7165, 2.7176, 2.7187, 2.7198, 2.7209, 2.722, 2.7231, 2.7242, 2.7253, 2.7264, 2.7275, 2.7286, 2.7297, 2.7308, 2.7319, 2.733, 2.7341, 2.7352, 2.7363, 2.7374, 2.7385, 2.7396, 2.7407, 2.7418, 2.7429, 2.744, 2.7451, 2.7462, 2.7473, 2.7484, 2.7495, 2.7506, 2.7517, 2.7528, 2.7539, 2.755, 2.7561, 2.7572, 2.7583, 2.7594, 2.7605, 2.7616, 2.7627, 2.7638, 2.7649, 2.766, 2.7671, 2.7682, 2.7693, 2.7704],
    s::Vector{Float64} = [0.2822, 0.593, 0.0691, 0.1474, 0.3586, 0.1121, 0.2773, 0.1418, 0.6625, 0.6, 0.4469, 0.8197, 0.5376, 0.0234, 0.7193, 0.5876, 0.7904, 0.6975, 0.8254, 0.4388, 0.796, 0.2523, 0.6688, 0.3571, 0.0893, 0.1872, 0.7914, 0.9879, 0.4024, 0.7375, 0.1406, 0.1612, 0.6001, 0.3847, 0.9423, 0.3092, 0.3507, 0.2995, 0.9904, 0.2004, 0.8696, 0.4609, 0.2579, 0.4172, 0.898, 0.7767, 0.6364, 0.7376, 0.2977, 0.8916, 0.8343, 0.9718, 0.3692, 0.973, 0.1912, 0.4646, 0.7714, 0.3854, 0.3021, 0.0944, 0.3479, 0.3945, 0.4981, 0.6592, 0.0782],
    t::Vector{Float64} = [2.73245, 2.75412, 2.72354, 2.70649, 2.70924, 2.7616, 2.72101, 2.70924, 2.70803],
)
    # 4.955496540068397
    T = m[2]-m[1]
    t_out = zeros(Float64,length(t))
    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i] += sinc((t[i]-m[n])/T)*s[n]
        end
    end

    return sum(t_out)
end

rozwiazanie()

##

          
function rozwiazanie(;
    m::Vector{Float64} = [-1.4, -1.3923, -1.3846, -1.3769, -1.3692, -1.3615, -1.3538, -1.3461, -1.3384, -1.3307, -1.323, -1.3153, -1.3076, -1.2999, -1.2922, -1.2845, -1.2768, -1.2691, -1.2614, -1.2537, -1.246, -1.2383, -1.2306, -1.2229, -1.2152, -1.2075, -1.1998, -1.1921, -1.1844, -1.1767, -1.169, -1.1613, -1.1536, -1.1459, -1.1382, -1.1305, -1.1228, -1.1151, -1.1074, -1.0997, -1.092, -1.0843, -1.0766, -1.0689, -1.0612, -1.0535, -1.0458, -1.0381, -1.0304, -1.0227, -1.015, -1.0073, -0.9996, -0.9919, -0.9842, -0.9765, -0.9688, -0.9611, -0.9534, -0.9457, -0.938, -0.9303, -0.9226, -0.9149, -0.9072, -0.8995, -0.8918, -0.8841, -0.8764, -0.8687, -0.861],
    s::Vector{Float64} = [0.6193, 0.5963, 0.4465, 0.6654, 0.2612, 0.3973, 0.6227, 0.8306, 0.8232, 0.1645, 0.1166, 0.1533, 0.1368, 0.3115, 0.3782, 0.5853, 0.9916, 0.7203, 0.4577, 0.9805, 0.5924, 0.3446, 0.3225, 0.5733, 0.4255, 0.0748, 0.3218, 0.3842, 0.9349, 0.3721, 0.4803, 0.523, 0.8254, 0.9765, 0.5556, 0.1769, 0.7406, 0.8094, 0.9228, 0.3776, 0.4377, 0.2669, 0.0434, 0.1523, 0.2589, 0.1039, 0.1144, 0.5643, 0.1295, 0.3887, 0.9947, 0.4753, 0.9483, 0.0581, 0.3499, 0.5155, 0.4584, 0.5849, 0.2568, 0.8713, 0.7422, 0.7046, 0.9644, 0.846, 0.7765, 0.3688, 0.7728, 0.5231, 0.9947, 0.7432, 0.0718],
    t::Vector{Float64} = [-0.93954, -1.15283, -1.09277, -0.99344, -1.15591, -1.2614, -1.12126, -0.90181],
)
    # 4.734115748228453
    T = m[2]-m[1]
    t_out = zeros(Float64,length(t))
    for i in 1:length(t)
        for n in 1:length(s)
            t_out[i]+=sinc((t[i]-m[n])/T)*s[n]
        end
    end
    return sum(t_out)
end

rozwiazanie()

        
        

        

        