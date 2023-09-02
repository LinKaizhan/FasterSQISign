#include "precomputed.h"

// each basis entry is a triple of the form P,Q,P+Q
// this is initializing the point using the classical representation {0,...,p-1} for elements in GF(p).
// We don't use this representation for actual computation but rather the montgomery representation (the conversion is made in init_precomputations using the fp_enc function)
// hence the ***_uintbig[] defined below should not be used in any actual piece of code.

const uintbig torsion_basis_sum_uintbig[3][2][2] = 
{ { { { 0x7d76bd8aa5b4284ULL, 0x13ae067bbcbfb67fULL, 0xf8167ff99caab13dULL, 0xf9553f385508903ULL },
  { 0x1b993d8fcf6153e0ULL, 0x93e09b3a94a7a063ULL, 0xb4ade22129e3d14aULL, 0x260a2dd7ecd2a4b7ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x40efa7f903343832ULL, 0xcd29fcaae15c2ee8ULL, 0x355654637eb313e8ULL, 0x32a2566712da0f85ULL },
  { 0x6f4f175e1b2dd70fULL, 0x4e243765c80da2adULL, 0x7de830a6ae361ea0ULL, 0x234d586d01f58b0eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe46948d83907e8f0ULL, 0x49df5dc9aa837112ULL, 0x3c608d3d1e3b0a2aULL, 0x2a7ba640f1939526ULL },
  { 0x7f92bad91d07988cULL, 0x12f76aa059a73745ULL, 0xc8358bc5702ab625ULL, 0x2f9c73b0ef867676ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };
const uintbig torsion_basis_twist_sum_uintbig[3][2][2] = 
{ { { { 0x1e799d85d670230fULL, 0x92bb5fd264da97e4ULL, 0xd76ec4948934c2dbULL, 0x147976c7d072afa5ULL },
  { 0x26430b79c0e1c0c0ULL, 0x1030d3ec57561d6ULL, 0xf5551e04f8b34951ULL, 0x1d22c5ec48c10d37ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa30b96f43988985dULL, 0xa63387d4da999955ULL, 0x910c7bc44c2029c1ULL, 0x166d9c5264f4100ULL },
  { 0xc84497f88517adedULL, 0x4bcbef9912256711ULL, 0x8eec62527dab25dfULL, 0x115435cb1983730fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7ef5663e6d8354b9ULL, 0x92727ff6e162d0bdULL, 0x905e99210401bff7ULL, 0x13c0fa0279eb1d32ULL },
  { 0x39a3832d019059ffULL, 0xbf4c9ec3e883869dULL, 0x662ef790edca435bULL, 0x29eeef656364966bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };

const uintbig torsion_basis_ted_sum_uintbig[3][4][2] = 
{ { { { 0x7a3f91bb49fe3b47ULL, 0xf7db234703509097ULL, 0x9b80a16ba1440e5dULL, 0x218cfdb71376f28cULL },
  { 0x3094e1ca7b6b2272ULL, 0xfbbed128b2d14307ULL, 0x3fc4e722163c0e4ULL, 0x30c3b8d8396a9905ULL } },
  { { 0xa5ec2dad745062b5ULL, 0x6970aba3d26d16a3ULL, 0x904c0a4e42d2f1aULL, 0x15199cd25dbc3cd4ULL },
  { 0x9906c014a20753b3ULL, 0x289d1d4b7ee35475ULL, 0xb02555a5ceb93326ULL, 0x2d64db808d5e99cULL } },
  { { 0x1ede022b32eca473ULL, 0xd159bda9a527296cULL, 0x840dc667253af0ecULL, 0x46bc58b65c35a12ULL },
  { 0x2804c51a49eb3446ULL, 0x8091a58e33b7b70eULL, 0xd02fb85eb051309fULL, 0x96e55c47985e8dfULL } },
  { { 0xb230b2fc36a956feULL, 0x3d42d41d83908d73ULL, 0x4b995790f01bb72cULL, 0x29412b5881ce9421ULL },
  { 0xaa1116a164fa8ac3ULL, 0xbf7e291104fe50ccULL, 0x4aaa313db11751fdULL, 0x15aa0f890eb84f06ULL } } },
 { { { 0x29ef58e53ca0373aULL, 0xcbf425c39d1339b0ULL, 0xb43c49200e8f258fULL, 0xd252a6d0ab300f0ULL },
  { 0x9e5b3be33fffb110ULL, 0x474c2fd78730e5b4ULL, 0x98af0ede9db0ef3dULL, 0x1322172cca2d224ULL } },
  { { 0x3701ea6bd7717338ULL, 0xdf6f84e4a778a771ULL, 0xa98383d57579462aULL, 0x50b30617769dc97ULL },
  { 0xd1abf3a4fad996a3ULL, 0xf0e6e784d1a1edfcULL, 0x1bf64dd2ae19dd44ULL, 0x9f836ad1e2d76a3ULL } },
  { { 0x834aa030aff7d650ULL, 0xbd7751fa23a424c1ULL, 0x3f3e1f428011f560ULL, 0x1973e6548ac7f18eULL },
  { 0x3a5e389c2e64f298ULL, 0x71842ed31223c404ULL, 0x837cbdbae5cfbc07ULL, 0xa67087d2d52d4f6ULL } },
  { { 0xc0eeef54050bd655ULL, 0x45dfb14918b96047ULL, 0xc743f16fa2944afaULL, 0x10cb0c5c2c2fa100ULL },
  { 0xcab61de1f6cc3c0dULL, 0x1df5fe8f059734c6ULL, 0xacb6e1493487d01cULL, 0x273038085c2fab6dULL } } },
 { { { 0x527a60a57015a46bULL, 0x5f0bee49aab3713aULL, 0x592fffb7166630deULL, 0xec94d2e914e0277ULL },
  { 0x709b62fbc66175e1ULL, 0x3e328b8c5b6c8e3bULL, 0x88d581a375ad2483ULL, 0x1b0249d2108b4ccaULL } },
  { { 0x705739ad18643b76ULL, 0xf47626d03e99068dULL, 0xe876bf786b301304ULL, 0x23b14f147191bb5bULL },
  { 0xc47d1c8462d46c65ULL, 0xa22236015cb5a966ULL, 0x27918f1114dc7811ULL, 0x321f0664ca13954aULL } },
  { { 0xf2078950621548bbULL, 0x978c620bee75579cULL, 0x9cc31ddf5e6d2899ULL, 0x9a5ad1fa6d87221ULL },
  { 0xc27247dbcf8e6543ULL, 0xbd8a31f793adb7e3ULL, 0x660fdffcf03062e7ULL, 0x1e18b149b7afc024ULL } },
  { { 0x3df2f6396df41b4fULL, 0xc5ca254d1f1165feULL, 0x25c631ae3467e8a5ULL, 0x2d4e2ae95f650f9bULL },
  { 0x52a957b4c09fc39bULL, 0xa5457567ce0c698ULL, 0x28618cb3889f8770ULL, 0x1ba85cad5b406505ULL } } } };
const uintbig torsion_basis_twist_ted_sum_uintbig[3][4][2] = 
{ { { { 0x29f0552b68aeb40cULL, 0x377df0b4aa4f3290ULL, 0xbe4ef0de6985cacdULL, 0x20f706975b2b98cfULL },
  { 0xf1801914bd40983dULL, 0x9eec9d5c3e29cb73ULL, 0x372f23bcfe5b5489ULL, 0x9a27471c43ff3e7ULL } },
  { { 0xa84548a492384b59ULL, 0x6c66ac7d6fbcbe2fULL, 0x389408cf4441ec1dULL, 0x5b780c8bc992a5ULL },
  { 0xeaa966416afb285cULL, 0xa954a6611832da99ULL, 0x787fce610f6b58a5ULL, 0x31c8da9230817c57ULL } },
  { { 0x80608c07c3011308ULL, 0xc4b5240014c2cf3bULL, 0x2df95db57b0b5244ULL, 0x11310f7739c50c2dULL },
  { 0xac0af12bbb86af30ULL, 0x19141767512118c3ULL, 0x28a0fdf13f7d7d61ULL, 0x18f3efb55266ab54ULL } },
  { { 0x884f3ba7d4d2fe8cULL, 0x6ee424ceb97b533bULL, 0x181277bc60550669ULL, 0x128321bec22eed2bULL },
  { 0x7c4d05cefbb2697fULL, 0xb662abd046cfed14ULL, 0xf731499ef02e52e2ULL, 0x81cc7ed39f877e0ULL } } },
 { { { 0x99498c3c902850f0ULL, 0xf2fd4cb1732c5260ULL, 0x7fd54860764a2eb4ULL, 0x1fef20e59a381a48ULL },
  { 0x4594e0664c99f434ULL, 0x941a9e786796775eULL, 0x2ad24cb5cb489dc9ULL, 0x288506c42dca8990ULL } },
  { { 0xcd245e23c2e89efcULL, 0xae7695c1872df118ULL, 0x3946ad919f9a6140ULL, 0x101b75dc57f8e47aULL },
  { 0x1b003725c3be026ULL, 0x6e190993fd6a008aULL, 0x1f94a20ed3df12dfULL, 0x281f9bf9b6e07c89ULL } },
  { { 0x62828cf4de9cfd2ULL, 0xc1d57c725050d957ULL, 0xa27941e9ec0f9accULL, 0x19db98875b72b7f7ULL },
  { 0x58d569177e7f2a2bULL, 0xec6571ab82d7db76ULL, 0xf94e969139db7e32ULL, 0xfdf1578f40c490bULL } },
  { { 0x9a242b4b0519a195ULL, 0x4f5fd2fd14722538ULL, 0x66fec4cb7c7d95f7ULL, 0x331bc70b5dec403dULL },
  { 0x2a8f81e8888eb6c9ULL, 0x9446c768eaf0bb3ULL, 0xad51ef220015fbc8ULL, 0x13f8801b4837b675ULL } } },
 { { { 0x3b6c218881850694ULL, 0x4f74fd2460d3c8aeULL, 0xdda732573bc65d67ULL, 0x7a3aca62464298bULL },
  { 0xb05e81a8ebb5eb86ULL, 0x5038b9aa7bef9c77ULL, 0x5316d050cb0f72b2ULL, 0x239e2e67a046979dULL } },
  { { 0xcd94561605174beULL, 0x2953801ea4305881ULL, 0x7a97f95a45e067b2ULL, 0x34356e134bd74b79ULL },
  { 0x564018d1165bf7eULL, 0x4a867151d4445e00ULL, 0x4a6ac519ab19fa47ULL, 0x3415af8ed0eb96a1ULL } },
  { { 0xd947e5c99a6fb0f6ULL, 0xae69bc7bcda703b4ULL, 0xc07f9d0612a39acbULL, 0x1bb7d68681f431f9ULL },
  { 0x59ea88fdeff96c0fULL, 0x74ece2b713c76091ULL, 0x86027e1e1450af9fULL, 0xd531dfb510842f9ULL } },
  { { 0x4f36668fe73e8b37ULL, 0x1fa64bdca1f9984bULL, 0xa1feaa85b07f552cULL, 0x18cc3b144695c746ULL },
  { 0x65a452b43ef634e6ULL, 0xc95b65778341075aULL, 0x89c66897021a207eULL, 0x3012013c32b1c06bULL } } } };

const uintbig torsion_basis_uintbig[13][3][2][2] = {
{ { { { 0x3afcf89589c297adULL, 0x2ad781b42dc1e57eULL, 0x34779003d35cc66bULL, 0x2aab2e54020be614ULL },
  { 0x4f24d3b7e79b6514ULL, 0x8377e1ae17725f41ULL, 0x8ada10ea1e69b751ULL, 0x107b4b2c2a7066fcULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xdf4ef839ae4ce3dbULL, 0xa240d353d3813e99ULL, 0x1692ce0c01e25861ULL, 0x3181b59dfb44a88aULL },
  { 0x623e5a391d5323aULL, 0xfb12f99fc94c1251ULL, 0xd873447945ac50d0ULL, 0x32efc25b40b1cbe4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x84087c8e0e4f8b56ULL, 0x683c4a41aa110fa1ULL, 0xf95f43eef87e9da9ULL, 0x215643acf8e38c04ULL },
  { 0x2e1b99b7ceade75cULL, 0xd021f63a7793f1e8ULL, 0xf3ef43042744c718ULL, 0x20137dcf80de1347ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x864cc9fbfa5d258ULL, 0x21634291fd177e32ULL, 0x2ffbb278d413055dULL, 0x1054ac2992e1e642ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xae96eac9f758e4feULL, 0xf7070b410c4840e4ULL, 0x2255699aa1b8c4cfULL, 0x2d4f53eaba693adbULL },
  { 0x5e0aa07b02ff4fd4ULL, 0x6a31cc2a46fd9c98ULL, 0x36739cd988951837ULL, 0x4c3252cc7360d07ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x34ecee36b7fe21d7ULL, 0x8b67a9761d4f9734ULL, 0xc557fec2c94ab47ULL, 0x199d99064ff365d8ULL },
  { 0xc1bf93f8f295d6eaULL, 0x5e4a6485cd08ff9bULL, 0x9349847a02a43ab3ULL, 0x2c1881a520261837ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xc4ea9138104c9bc3ULL, 0x77d0f93925a1722ULL, 0xf585c4e6cb79f788ULL, 0x2aedb4427071535ULL },
  { 0xa1f226d8d981a148ULL, 0x629f4e6766ab6908ULL, 0x8b112631fa0bec77ULL, 0x307515607b094cdaULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x116f33482125fd51ULL, 0x59d2aaa447ebb062ULL, 0xfe3be99d0e3733fdULL, 0x291ad2a7d1670c3aULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0xa3570bd8f858acfULL, 0xb828e7ef8484d97cULL, 0xc4f4367d6ffdee95ULL, 0x1cbfe3e0ac467ea8ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x3799fdc3ff2249e3ULL, 0x7fcacf89761ccdc1ULL, 0xb5c7d96805e79fa8ULL, 0xce0052a7699bc3eULL },
  { 0x55dbef22531525dbULL, 0x557574e178e624dcULL, 0x6cdaecc77312a345ULL, 0x310b7a17d67d10d1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x21562c24a8fd73c4ULL, 0x3fe5eac3b161e495ULL, 0xfb53b5ce6a46f4a4ULL, 0x6b2e86b40ccd85eULL },
  { 0x3c2afb8b7eafcc9cULL, 0xdd11ef47bd2ed2bdULL, 0xc420e00cce396e61ULL, 0x1ad6839210d06323ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x829b9d37bcb6d348ULL, 0xbdf2673cc7aa8affULL, 0xf4c630c5d9272db5ULL, 0x22dd47af93d04639ULL },
  { 0xeca6d1eb7f5bfc07ULL, 0xe6f657df05cd5a42ULL, 0x5838549853703d08ULL, 0x870a4232948d5d8ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xe72fb4edb64cb143ULL, 0x3d94a50ca4c97d45ULL, 0xd8467ded31c79126ULL, 0x31693f9ab65e77bbULL },
  { 0x83b6a361046c0bb5ULL, 0xe4784483c03dc4feULL, 0x7d33695f6f12d116ULL, 0x2d3f4458183e0610ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xaabc3d21b640dcdeULL, 0x781c841011eebe04ULL, 0xe99f4f6f8f340128ULL, 0x24b028b3c016d1b7ULL },
  { 0xc4302886028a13e8ULL, 0x2ee4f36217e58fc3ULL, 0x2963b9d6b7c375cULL, 0x256c8a2b9019dc67ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb58b502688945d8fULL, 0x533d6a2f2803df81ULL, 0xe0f9e75b924ee0feULL, 0x31bc2680300102bdULL },
  { 0xefab6f00eae38ab1ULL, 0x6956885c6c5c3873ULL, 0xb5341203e2112f8eULL, 0x20e5cdb0b0aacc0aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x381ffd4929dd20bULL, 0x8e8a1353337e02bfULL, 0x6771daa9e4495ad6ULL, 0x6adbdde91c3676eULL },
  { 0x72ddf7eb8622e0c7ULL, 0x6ad30d40065df54cULL, 0x39d91aee9e7e9537ULL, 0x2a6a02b51d825bdeULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x983dfcd09fec32b0ULL, 0x8176315cb26a1b5bULL, 0x6ab036f05692a32eULL, 0x2a04e88a20ab0fc6ULL },
  { 0x48855e4e9a705871ULL, 0x64eaeb3b3fa2f8dfULL, 0x6546fef3fd6348a2ULL, 0x306d039d59ba0c21ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9701ab915efd2c93ULL, 0x5596fb5505ed15f6ULL, 0xc0b14dd9553997efULL, 0x1c00eae5c3df7570ULL },
  { 0x4a6119f260269e5cULL, 0x8f09039a7d6e79a9ULL, 0x71ed394c8ef58bdbULL, 0x9d2d533fa505542ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xc200c3d3080a2448ULL, 0xe0471425cd3f4be5ULL, 0xdcd2bc68c0df64dcULL, 0xf0888af6f9c32fdULL },
  { 0xb510b713fd598984ULL, 0x41f46aff91edc747ULL, 0x9aad893f8bfd4e9bULL, 0x33f7717878c13e77ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa31f09b9f98e508fULL, 0xb670e71a2a5d5645ULL, 0xf788db3a394c318cULL, 0x95fb93fc1a02955ULL },
  { 0x318a9b1522663040ULL, 0x77949b1a72a27236ULL, 0xddeb6807f45d3248ULL, 0xbd74d1b1ac09d9bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xeff073cdce2fd579ULL, 0xf8b1869e9b0a7520ULL, 0xcf27d49c97c0395ULL, 0x42b08b10c0816a0ULL },
  { 0xa73fd70a6a1c36cbULL, 0x61ee56fed5c27b3aULL, 0x257731a7f771ce73ULL, 0x1115c3167e02cfd7ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x98c813d03dc51f09ULL, 0xada8609e199acca4ULL, 0x364d4a66c2bbbcb3ULL, 0x106a83a72f362707ULL },
  { 0x1d3f2c6d3374e833ULL, 0xebada9df33dbd6baULL, 0xe574eff467678db5ULL, 0xcd181668532eae1ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x79f5bfb176fdbd3dULL, 0xa4b7c87630c29205ULL, 0xfa64e4add7fa3ae1ULL, 0x1f0d5a944d836b8fULL },
  { 0x6026ee2858d14345ULL, 0xf94d697f636fcac8ULL, 0x296f5064eaba2c08ULL, 0x15e96819c6dad279ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x53c71cfe07ec6732ULL, 0xf464786d4c91041dULL, 0x69a65b96a110781eULL, 0x2d7172dee7172a32ULL },
  { 0x7a64c8ae74460f51ULL, 0x96887fab567fa9cdULL, 0x6e20207267382296ULL, 0x7639fd8ec825ff3ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xdffb36fe8a7e1aa4ULL, 0x50f7737f7a317496ULL, 0x2c5a5367e718193eULL, 0x6094f45d52ee758ULL },
  { 0x787a77fbe3d06992ULL, 0xcd15241f5b99c714ULL, 0xe6bef985bb2282aeULL, 0x115773a30adbdea6ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa2940c737d0f1f0fULL, 0x434e45195ea4012dULL, 0x5d043b0bf279971fULL, 0x10ff98de4b4b7b7cULL },
  { 0x5b8aad1ae08766aULL, 0xd239922e01bafb6ULL, 0x4ec5d70a00e5864eULL, 0x2a68c0114080aa01ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x938d28315817030dULL, 0xc8a6b156f6a2c74aULL, 0x5f2421e3e091ea7aULL, 0x1a70315f874cc2fULL },
  { 0x6bcdae3acc3c56feULL, 0x7cdf4fbb6a400ccaULL, 0xdf9f741947b92986ULL, 0x213227233e111556ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xe7e7f868850055aeULL, 0x4d85b4aa2b87482aULL, 0x1106512703441ae4ULL, 0x22736304c85e462cULL },
  { 0xa3cdbda4ba868493ULL, 0x1f0320337b2e6d8dULL, 0xd9d9bf666c084eb4ULL, 0x1e94d908177fe9aaULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x981e9757df61d8ecULL, 0xc2ad565453a72e06ULL, 0x268358208a71abe5ULL, 0x188801e8aed3eabfULL },
  { 0xb2a6a08551f20605ULL, 0x304b8219b38bfbdaULL, 0xe7223f3e291289c5ULL, 0x1ad8ac86ddfa272eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xeddc3fa8c2377698ULL, 0x1a7fbb4eb830640aULL, 0xeb731b436613a090ULL, 0x301a89792346e3beULL },
  { 0x11fc13dd984d3cd1ULL, 0x1fa58db2155dd62dULL, 0x95c5e662454488bfULL, 0xf4698c2a79256a9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xd2c482e08b26cd4bULL, 0xb66355de67ab47eeULL, 0xf120975466f624a9ULL, 0xddb3ee068fbc7afULL },
  { 0xcea66431a3eb1c21ULL, 0x1beaed6cc4c54e27ULL, 0x6caa6c229f50d8e9ULL, 0x626d3b74db836bdULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6846defda3cc1d76ULL, 0x65ca969b43d6b27eULL, 0x32a38356fb4f0949ULL, 0x4e64065a7071054ULL },
  { 0xccc06d165c497daeULL, 0x594eed3ff3a0367cULL, 0x74d49e398e1d58acULL, 0x4fd3358c013ed5aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb789126ad6ed0e82ULL, 0x3c69dfbb226661f9ULL, 0xdd533e6b358dc7c9ULL, 0x2cab3427d0384ae7ULL },
  { 0x1c0d843bf20893f7ULL, 0x4b05110ca0b85583ULL, 0xa789455d72d8a8f5ULL, 0x285f545415780401ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xdaf1db136d88a344ULL, 0x4a80e3919b0492c8ULL, 0x9277cdcbe7080df7ULL, 0xcb2d84da9203077ULL },
  { 0x2f4280246e334984ULL, 0xde24759a5ea7545aULL, 0xb9b2d65ac8c2d356ULL, 0x22edb209eb6f4ed3ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x36fc01abfbb42c69ULL, 0x4ed9dd6dfd108c5cULL, 0x9a8ca3b070ecab46ULL, 0x2661cc99e525bc9eULL },
  { 0xf085a62fdbef5df8ULL, 0xf76d43df34291bbULL, 0x5811a7b3ca09ca1cULL, 0x1db3f947c57f5e6eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x616b2428cac07eb9ULL, 0x6484bf94300354ceULL, 0x47bf9617b9092c97ULL, 0x149bc5332eb89e61ULL },
  { 0x5fe51bc9db8c3af6ULL, 0xcd4d15253d0003fcULL, 0x8325e9d47fa3859dULL, 0x155cde23d5acfc27ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x1c222162274d2852ULL, 0x614b6b2fb92fbfdaULL, 0x4e3251b2e14b6e20ULL, 0x315ffee96682da49ULL },
  { 0x99255684d92053c9ULL, 0xb21a90d5ce05ff5ULL, 0x2c46f1896db7ea0dULL, 0xec2b9a1e081eb98ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xba21bab7a83224a6ULL, 0x3fe104f8cffce7fdULL, 0x6ee43f970aa54ca0ULL, 0x344a19cefa3d3963ULL },
  { 0x8c3327191317bf36ULL, 0xaf7b36b6f355af9dULL, 0xbd6b9ab584e7f1c5ULL, 0x9377f02ebff2575ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x42a65e32b4d6b9c7ULL, 0x3690e8f2d0df0c81ULL, 0xac6823551fdd3b31ULL, 0x1f39473a02551ca3ULL },
  { 0xc3be21e261591248ULL, 0xff52c9d6d28500f2ULL, 0xfb573012018d4080ULL, 0x3284f67862662d0eULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
};

const uintbig torsion_basis_twist_uintbig[18][3][2][2] = {
{ { { { 0x8a611fb64be5e963ULL, 0x13877c1a6164b6d0ULL, 0xdf5878a7c50cc63ULL, 0x18f980ff9b7d78b4ULL },
  { 0xc71d55720d69c72dULL, 0x19070c2ff840fc27ULL, 0x5b5d7691b4493b0ULL, 0x7196756d36460bbULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x83d9e9dbef81a5a8ULL, 0x48f912991d959382ULL, 0xf18acfc34f4f3209ULL, 0x8acf369f133d511ULL },
  { 0x7e0b6d9db0a511ebULL, 0x34996ded7d23038eULL, 0xfb03323602e137f6ULL, 0x14419a6a72588f0dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6bbefc5069781544ULL, 0xe9d3527e63fc55f1ULL, 0x54ff61cd9a4618feULL, 0x1a5fd516f68e98f0ULL },
  { 0x8e4ef2516e612a44ULL, 0xedc3100806dab231ULL, 0x9daf2157e56c1e2aULL, 0x2b6cc445867d3ba6ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x7714ff911f4c0314ULL, 0x2558921b59e41a51ULL, 0xb7444c778fd8e3bcULL, 0x2683fdfd56c8256bULL },
  { 0x87acb7b90f935d7fULL, 0xd65571177bc53b9fULL, 0xc60726b3fee79111ULL, 0x21a29ebb6e4528f9ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe2a65d86743be51fULL, 0xacf73de94f44efabULL, 0xab7ca669405c1f5eULL, 0x213f4dd0241ee620ULL },
  { 0x7e213ec2ddb7e1aeULL, 0x8b219323c3bc3a8ULL, 0x7a7e43334cad21c7ULL, 0x1e5a675579003dacULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x65c2a00f4857e799ULL, 0xb1d21f1135d8b6eULL, 0xc8aa5f27933b2ed9ULL, 0x13d18577aacf94dfULL },
  { 0x8dd1a1cd39e2737fULL, 0x588aa2adbfc0ab91ULL, 0xed123901877139e3ULL, 0x173c842ba3aa391ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x934cf14f44713b8cULL, 0x64448bdf15cf22fULL, 0x2e36535408dcf2dfULL, 0x1dd5a14f3890273eULL },
  { 0x822a97cf28d8b85bULL, 0x69cb63138c88734bULL, 0xeb0d19ae3bdf1c83ULL, 0x28038cf82929837aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xaaca4ed1036d28fULL, 0xbd9ae954af734a98ULL, 0x6f87d5a8a321b273ULL, 0x1d2e9a510f099b03ULL },
  { 0xe9c668f908e0f7c5ULL, 0x7252087ce7a550f1ULL, 0x9811bb9480a53e27ULL, 0x163a62d9bc07f222ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x92ca7aa2ac46cca4ULL, 0x14b461503a2437cbULL, 0x8a2af35af639d05fULL, 0x30627d7973160ce7ULL },
  { 0x7b466856cb061378ULL, 0xb774f5b15c60a9cULL, 0xde6d37895d119326ULL, 0x17045ff973137182ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xd728e611780d685fULL, 0xb27882ea6bcb6e33ULL, 0xfabcd6f36c3d6e18ULL, 0x4f6d3816218731fULL },
  { 0x2863b98ad5a70b07ULL, 0x4ef6a7ec863ffe98ULL, 0x8d496e225384a51eULL, 0x27ca042905f4ba00ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x84d0b51aad28aaceULL, 0x3400185075b520e9ULL, 0x9379ad66af4d16d2ULL, 0x2b04f1b452d29b9fULL },
  { 0xd4755e8f2ea33b88ULL, 0x51f483ef43b36a6ULL, 0xb60bb78bb52f551fULL, 0x3142d38013a8a248ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x577f2540c95bd9d9ULL, 0xc56e5d610ae0034bULL, 0x35bef2f4718f3ed8ULL, 0x130db5bc4aa4cf43ULL },
  { 0xcfd689b84cbaf664ULL, 0xda43765ccb2299dcULL, 0x5c1d962ca63a2469ULL, 0x16b2b1049f5c5a2cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x98d04ea7d609412eULL, 0x6a1ee05e867c7a89ULL, 0xa2f7c9550332bea8ULL, 0x3286b72473818033ULL },
  { 0x8482f8272c29e462ULL, 0x71f08f17137d8d28ULL, 0x3c6b4e8a8a6432acULL, 0x276ba56ea2d3972aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xe9ae927e2ed8a428ULL, 0x2905a548f2b975d3ULL, 0x5bc06cd8bf6c2d16ULL, 0x1cc1084369453c39ULL },
  { 0x1e04603ccd5560d7ULL, 0x5a0dadf3297fdcc7ULL, 0x368091a8267b2edfULL, 0x24cc9a8857f36a8fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x3c1579112f4b0e59ULL, 0xf38be7fda4d9216dULL, 0xe31b495535f69f1ULL, 0x27d9c56a85ca0044ULL },
  { 0x1f270afc7837990dULL, 0x20451881f2bb8c25ULL, 0xe6662877ba80269eULL, 0x26c41cab2593dbfeULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xd55362b3f39dca32ULL, 0x7ae75adfff7a64d9ULL, 0x9f7d65cf2826e0a9ULL, 0x238b669fb1932544ULL },
  { 0xc34807d166049a46ULL, 0xc89ef732586d30faULL, 0x5b9cc3caf168e59fULL, 0x3ca72bc8d32c266ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x7f59542dc0fbe112ULL, 0xb82f0a4857555d3aULL, 0x5fa23a5782f233d8ULL, 0x2153596b9b22c5caULL },
  { 0xa0b0951339a6d74fULL, 0x20995eea69163babULL, 0xd5d5d4e088ffdab7ULL, 0x2fcf5236856b085cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9e332a8b74aebcd7ULL, 0x55ee7dca776981d8ULL, 0x9a323b325406feafULL, 0x373f6d5b8666bb1ULL },
  { 0x7d6f157c28bea547ULL, 0x3d12bc00cb52ddc4ULL, 0xa0879c9ec06b6c81ULL, 0x1e2ba7d054c12620ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xa3fbd88ae14a27e4ULL, 0xb52b783b8ba37072ULL, 0x7b85837be562c7c4ULL, 0x313006f2f31f129ULL },
  { 0x4fdd7e7c963575dfULL, 0xd7e813ba1ca47460ULL, 0xb3d1a370a6d87c4dULL, 0x11eb133a98165f1aULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x9eccf8464244f4b5ULL, 0x3db33671cb42c2eeULL, 0xdd16cdff0d5250cbULL, 0x16947b08680c9e5fULL },
  { 0x3d1884563bac43c6ULL, 0xe37c8f240f0fb0b9ULL, 0x3cecf16713f445c8ULL, 0x88cfb19d425ee0fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xb658f56ad94ecaf8ULL, 0x538bf3b25837afa4ULL, 0xbbba63622f2764e0ULL, 0x23f825283fbc8f07ULL },
  { 0x613a0b3b12265983ULL, 0x8776c13baaecd098ULL, 0x1ec98cd32eba5f15ULL, 0x1baba4c848949c3fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x2427b99e2e2fd6baULL, 0xf94282ac382c70d6ULL, 0x5ec561bff56c4ad0ULL, 0x7d5078b401527b0ULL },
  { 0xfdfcaa4c241c0e7fULL, 0x502ed7808d715e77ULL, 0x51fb70f4db034e60ULL, 0x1438051620958e45ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x17ab4f2f09bd2ec6ULL, 0xae76480a0d904611ULL, 0xed463c3cf129662bULL, 0x602c25dd4d677d5ULL },
  { 0xe2c6cecd32f877f4ULL, 0x42929c9d5fcfa15dULL, 0x8cbafbfc79da4096ULL, 0x2c09f6639931a32cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xabcf83d43a3bbaf6ULL, 0xa8dd026e1e7e7350ULL, 0x6dd7069832a8c415ULL, 0x2116bbaae7694796ULL },
  { 0x746c5ea8d5d64ca1ULL, 0xe02a3f881d5de8c6ULL, 0x50a562c2b19125ffULL, 0x481d01c01b4caefULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x48fa481a09d97303ULL, 0xbc782862f01040a6ULL, 0x1ac55e84eb51b6a0ULL, 0x2563e64e230d41f4ULL },
  { 0xad5e4d2d1b3a4a7bULL, 0x351f053e009b34daULL, 0x12aa963365a30da5ULL, 0x218a43986b4b5d3cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x75258ebde5208683ULL, 0x1177433e16a34393ULL, 0xcc5b2fc2cc417645ULL, 0x2f7e38011a471353ULL },
  { 0x2b606c6453837ba3ULL, 0xf09fc92828c4d1f1ULL, 0x72db2a359f659fe5ULL, 0x319983c18ea66380ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x452b2bfe5f197830ULL, 0x39ad2fe0b67d71b0ULL, 0x54ab01e39a26bbeaULL, 0x21c5ca98c11ef5a7ULL },
  { 0xd9e39540ebc72c84ULL, 0xb965bf2ce9f853b0ULL, 0xc5406774407afe3eULL, 0x20e7f14017f88f6bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xa78471aed6013d40ULL, 0x4bcce1dec2d30d62ULL, 0xc5ca55f2d25786c5ULL, 0x6378305e3e1a0b6ULL },
  { 0x3f7740c6469ee237ULL, 0xc5c5b27f8539c3e5ULL, 0xe3335e295cbaa475ULL, 0x2433006f516d1b66ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x5eb0c51f2760eda4ULL, 0xeda690d39fc051ddULL, 0x69776e9fa17de33fULL, 0x14025be419e77c77ULL },
  { 0xc8aff83b712255aULL, 0xc10b9a397fd0b0dbULL, 0xcd2c9666306cf97cULL, 0x1027365a5c35a5efULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x53c0457e2423e4dfULL, 0xafc74f2bae89866fULL, 0xa0ea10ecc5ff2f1dULL, 0x19b5e3197f56ca10ULL },
  { 0xd7018caa8a33bde9ULL, 0x20979c0338503b95ULL, 0xade4aa544d100cd7ULL, 0x327d4b234b238561ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x14a394a57df331cbULL, 0x4af1883059e4bc40ULL, 0x67e651a125f25749ULL, 0x209f89adbd120b7fULL },
  { 0x4a106c5800e81937ULL, 0x3b653426ef2070f6ULL, 0xab9522ecfa60e02fULL, 0x2ccfaad3bc299368ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6100adaba2d3f026ULL, 0x47e25a9b8f83eb2aULL, 0xb503434287589624ULL, 0xb8edd28b09610bULL },
  { 0x427e29c27d7a7ad2ULL, 0x801210554ac6856bULL, 0xcbc0654b073e03b2ULL, 0x29a540ddac24e71cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x297350ec7322749bULL, 0x46843898358efc72ULL, 0x2626714a4a4dd296ULL, 0x1ee538bd1846353ULL },
  { 0x5375f69515efa386ULL, 0xac6576972368d6c0ULL, 0x59d0556f14924f9aULL, 0xc7ea59919c87557ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xa038114cd0292699ULL, 0x516b0d8146e3ac3aULL, 0xbdfd70b0e2c8c315ULL, 0x2e649a878664d7f5ULL },
  { 0x8391450a6ac2b9f0ULL, 0x624f828b8e1529caULL, 0xaea0cc2ee5690ab9ULL, 0x24d65adbf2e25c18ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x505bbe1d727b662ULL, 0xd418cd0e425a30aULL, 0xa0532673e61fcedcULL, 0x1aa2700cb189e869ULL },
  { 0xecf42f6dc010b942ULL, 0xb2ca8d30b0842a56ULL, 0xd97b1fdea860f395ULL, 0x5c1ddd2a607bed4ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xa54064ee638ac6fdULL, 0x3d039facafca3079ULL, 0xc66681beddc1b666ULL, 0x292c89d9ed39cb1dULL },
  { 0x626c6610c989b7ceULL, 0x6414e379523d20bULL, 0x89ea55772ecf68afULL, 0x10a11a118bafa710ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xc7c5a5b61eb5c224ULL, 0x83c11e21e36ffc79ULL, 0x870359e341119357ULL, 0x5a9c31b272c79c7ULL },
  { 0x9af7836264663d44ULL, 0x2e0cb6978feb134fULL, 0xb9dd3b7c5fd0a701ULL, 0x5da7ef080576e3dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x6893f25725122da5ULL, 0x4ba0b21a280d783bULL, 0x90ebf13d0160d0f4ULL, 0x26a7c3356f852cfeULL },
  { 0x5e23f55bd861751bULL, 0x32cce04fc1063896ULL, 0xccab6a3ceac36154ULL, 0x287281a35a184673ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xf66277b9528d445aULL, 0x5887f6fd3c924924ULL, 0x2ca7178c180f72acULL, 0x1bdda4ef68bcaf5dULL },
  { 0x9567882772da747dULL, 0x4efae61d242de02fULL, 0x971d7a0d51837979ULL, 0x21e0b6f8560f7493ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x6e2c51cf0bf203c0ULL, 0xae8826db9e95c805ULL, 0x3faa0b0c7eb762c7ULL, 0x4b92d833051fd39ULL },
  { 0xe7a38b33d6a0b2f6ULL, 0xbb23e32071425d80ULL, 0x85a8e788a1d00ac8ULL, 0x158c42b326ae5858ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x66cb1876420eadc0ULL, 0x808aeec8089b9d76ULL, 0x51991c416dd7bd8dULL, 0x17901306b03dbf3dULL },
  { 0xaad9e84ed3394d69ULL, 0x9fa804aac8f0238bULL, 0x9184929b3ab6c61ULL, 0x27d716e34c8d4b1dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x52c36cc7248e2535ULL, 0x5840ef0ed1246bd2ULL, 0x7a727568ecd06229ULL, 0x1a0aaa0621328b9cULL },
  { 0xf556131f305115d7ULL, 0x64a0490f0451c593ULL, 0x8a5e559ddeef941dULL, 0x2cd66f8f432e9bc2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xb8752519dd64a1b9ULL, 0x14caf564c90d501bULL, 0x9539591437d1c020ULL, 0x1b160b8ab87634a3ULL },
  { 0x87293f58fa7a312bULL, 0x8205cd0de508ed10ULL, 0xa5244ee73a0e5423ULL, 0x1e74eb84a147318bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x8c841e1378e27ad8ULL, 0x3f1ea91d2898daafULL, 0xb09f1799a7efff4fULL, 0x304f969ff0313085ULL },
  { 0xa9d6643a8d7f2719ULL, 0x79126c952b8b90e2ULL, 0xeebe038e7c466481ULL, 0x64dd34824e0a618ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x56be47f12170e52cULL, 0x7e96059c8aa8d4aaULL, 0x43220db319e6328aULL, 0x33bd04a516e1b883ULL },
  { 0x2d77dc3125a35f5dULL, 0x40cd8df86fc6e16dULL, 0x5cea3c631ad43cf7ULL, 0xe827b1b57bdda2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0xff5d556764dab4b1ULL, 0xfec342ffdf0e52d1ULL, 0x72a29dea772dc64bULL, 0x3406ad9cd9355fa7ULL },
  { 0xb5397283ff711422ULL, 0x510ec74d21153e98ULL, 0xac81056e11bc517dULL, 0x129fd42b8613993bULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xfc227754a2d048e7ULL, 0x9818119342ad7298ULL, 0x77792b46c70f92dcULL, 0x2392c0985a532ad3ULL },
  { 0x2932a759c043358dULL, 0xe449d5c4c8e8867eULL, 0x21bb8e926ec170dcULL, 0x409061dc1bc9efaULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xebe62b8af44fcf94ULL, 0xfdb69e8bf4682c2ULL, 0xacb676e954e8c310ULL, 0x325c9770aea408feULL },
  { 0xb211722a2180eb3ULL, 0x2a48348dfb6d49e1ULL, 0xfb7ef7e955c87754ULL, 0xbad6b3c25e85e06ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x403f77a5e97b108bULL, 0xcac82fffcf463c14ULL, 0x9b305b932c560caeULL, 0xffe17c7b633eeeaULL },
  { 0x3d93d0f6d87795e2ULL, 0x74e3eae8ce779703ULL, 0x5c2ddd4f82be8627ULL, 0x4ce7b60db40bc5dULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xab44fdf1d31031cdULL, 0x811732328b6ff5b0ULL, 0xe073f78c646af83bULL, 0x2af3669f24914496ULL },
  { 0x8faa31ddbb784c94ULL, 0xe81ebec66fd2263cULL, 0x2f6777b9761617ebULL, 0x229f2ac0c0b39266ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x4cd5fcc235c18f95ULL, 0xdae556d95fc96c12ULL, 0x6a21842c4e8b7769ULL, 0x2fcdd68b79d93991ULL },
  { 0x3bd4b2957433513fULL, 0xd73240cc327d705ULL, 0x86bfe711675129fbULL, 0x2beaf4ef8f61a5c2ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
{ { { { 0x7b612b9cdb93c3e8ULL, 0xbfa8d398e5d7106dULL, 0xf6ddd02eb7ac540eULL, 0xbac4448f9c60a0ULL },
  { 0x413f232890d4ee32ULL, 0xa9754baae0e7b4beULL, 0x67ecf6d69f89e58ULL, 0x2ca8f93117bbe878ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x3f0537e1aae0013bULL, 0xb7943492016e5983ULL, 0x8d9787f325c535a1ULL, 0xd70e60c70347f53ULL },
  { 0xbe33eba02831d313ULL, 0xb16c1b347bb6c28eULL, 0x243aceaec3a3f140ULL, 0x1180afb0cdaf7b10ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0x37f12b71ce302af3ULL, 0x7f988816a1563102ULL, 0x6f7583aed0cb7cc0ULL, 0x32b67d36d513f756ULL },
  { 0x9216e6afcb19b6ccULL, 0xf1f951dbbac57441ULL, 0x20920c15d72799a7ULL, 0x22caa8186cbf0b55ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } },
};

const uintbig torsion_basis_two_uintbig[3][2][2] = 
{ { { { 0xb60a2ea27d0b92bcULL, 0x781d1698c4a7d3efULL, 0xb0f83c00e53f02ddULL, 0x9e39507ab7883d0ULL },
  { 0x1dae14ff48809b4eULL, 0x5af43e3929fc128ULL, 0xf6d180e31d3eaef0ULL, 0x31a276fc45de603fULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xdd1cfae93b4a2698ULL, 0xacb717de868ac38ULL, 0xdd91ff638ee5c4f8ULL, 0x2c04565ad62ba9c6ULL },
  { 0x96e9ddbb84d9d23eULL, 0x8db179a305b02fa4ULL, 0xbde759fd18122ac1ULL, 0x81daa6a2ddfd0e7ULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } },
 { { { 0xda1e9846a23ed1e7ULL, 0x21daa769ba502517ULL, 0x9d55a71b7d1cf80ULL, 0x1890feeaca8a8b83ULL },
  { 0xbb81dacfc2ce89e0ULL, 0xb04bff8e45721c43ULL, 0xbd2aa314162d32a1ULL, 0x30cdc51abf4f642cULL } },
  { { 0x1ULL, 0x0ULL, 0x0ULL, 0x0ULL },
  { 0x0ULL, 0x0ULL, 0x0ULL, 0x0ULL } } } };


proj torsion_basis[13][3];
proj torsion_basis_sum[3];
point torsion_basis_ted_sum[3];
proj torsion_basis_twist[18][3];
proj torsion_basis_twist_sum[3];
point torsion_basis_twist_ted_sum[3];
proj torsion_basis_two[3];


void init_precomputations_generated() {
	global_setup.action_2 = malloc(13*sizeof(GEN));
	global_setup.action_3 = malloc(13*sizeof(GEN));
	global_setup.action_4 = malloc(13*sizeof(GEN));
	global_setup.action_twist_2 = malloc(18*sizeof(GEN));
	global_setup.action_twist_3 = malloc(18*sizeof(GEN));
	global_setup.action_twist_4 = malloc(18*sizeof(GEN));

	global_setup.action_2[0] = mkmat2(mkcol2(stoi(22ULL),stoi(18ULL)),mkcol2(stoi(5ULL),stoi(3ULL)));
	global_setup.action_3[0] = mkmat2(mkcol2(stoi(11ULL),stoi(15ULL)),mkcol2(stoi(16ULL),stoi(15ULL)));
	global_setup.action_4[0] = mkmat2(mkcol2(stoi(5ULL),stoi(0ULL)),mkcol2(stoi(3ULL),stoi(20ULL)));
	global_setup.action_2[1] = mkmat2(mkcol2(stoi(4ULL),stoi(6ULL)),mkcol2(stoi(3ULL),stoi(3ULL)));
	global_setup.action_3[1] = mkmat2(mkcol2(stoi(0ULL),stoi(0ULL)),mkcol2(stoi(3ULL),stoi(1ULL)));
	global_setup.action_4[1] = mkmat2(mkcol2(stoi(4ULL),stoi(6ULL)),mkcol2(stoi(2ULL),stoi(3ULL)));
	global_setup.action_2[2] = mkmat2(mkcol2(stoi(7ULL),stoi(1ULL)),mkcol2(stoi(5ULL),stoi(4ULL)));
	global_setup.action_3[2] = mkmat2(mkcol2(stoi(4ULL),stoi(4ULL)),mkcol2(stoi(8ULL),stoi(8ULL)));
	global_setup.action_4[2] = mkmat2(mkcol2(stoi(3ULL),stoi(3ULL)),mkcol2(stoi(8ULL),stoi(8ULL)));
	global_setup.action_2[3] = mkmat2(mkcol2(stoi(15ULL),stoi(15ULL)),mkcol2(stoi(9ULL),stoi(4ULL)));
	global_setup.action_3[3] = mkmat2(mkcol2(stoi(16ULL),stoi(4ULL)),mkcol2(stoi(16ULL),stoi(4ULL)));
	global_setup.action_4[3] = mkmat2(mkcol2(stoi(5ULL),stoi(6ULL)),mkcol2(stoi(18ULL),stoi(14ULL)));
	global_setup.action_2[4] = mkmat2(mkcol2(stoi(708ULL),stoi(96ULL)),mkcol2(stoi(464ULL),stoi(133ULL)));
	global_setup.action_3[4] = mkmat2(mkcol2(stoi(439ULL),stoi(254ULL)),mkcol2(stoi(329ULL),stoi(403ULL)));
	global_setup.action_4[4] = mkmat2(mkcol2(stoi(109ULL),stoi(701ULL)),mkcol2(stoi(199ULL),stoi(732ULL)));
	global_setup.action_2[5] = mkmat2(mkcol2(stoi(806ULL),stoi(887ULL)),mkcol2(stoi(953ULL),stoi(563ULL)));
	global_setup.action_3[5] = mkmat2(mkcol2(stoi(1365ULL),stoi(178ULL)),mkcol2(stoi(246ULL),stoi(5ULL)));
	global_setup.action_4[5] = mkmat2(mkcol2(stoi(45ULL),stoi(51ULL)),mkcol2(stoi(524ULL),stoi(1324ULL)));
	global_setup.action_2[6] = mkmat2(mkcol2(stoi(6ULL),stoi(42ULL)),mkcol2(stoi(45ULL),stoi(41ULL)));
	global_setup.action_3[6] = mkmat2(mkcol2(stoi(32ULL),stoi(23ULL)),mkcol2(stoi(10ULL),stoi(16ULL)));
	global_setup.action_4[6] = mkmat2(mkcol2(stoi(1ULL),stoi(11ULL)),mkcol2(stoi(17ULL),stoi(46ULL)));
	global_setup.action_2[7] = mkmat2(mkcol2(stoi(127ULL),stoi(34ULL)),mkcol2(stoi(105ULL),stoi(70ULL)));
	global_setup.action_3[7] = mkmat2(mkcol2(stoi(51ULL),stoi(53ULL)),mkcol2(stoi(108ULL),stoi(147ULL)));
	global_setup.action_4[7] = mkmat2(mkcol2(stoi(102ULL),stoi(106ULL)),mkcol2(stoi(110ULL),stoi(95ULL)));
	global_setup.action_2[8] = mkmat2(mkcol2(stoi(175ULL),stoi(43ULL)),mkcol2(stoi(89ULL),stoi(88ULL)));
	global_setup.action_3[8] = mkmat2(mkcol2(stoi(171ULL),stoi(144ULL)),mkcol2(stoi(94ULL),stoi(93ULL)));
	global_setup.action_4[8] = mkmat2(mkcol2(stoi(40ULL),stoi(6ULL)),mkcol2(stoi(84ULL),stoi(223ULL)));
	global_setup.action_2[9] = mkmat2(mkcol2(stoi(19ULL),stoi(271ULL)),mkcol2(stoi(261ULL),stoi(262ULL)));
	global_setup.action_3[9] = mkmat2(mkcol2(stoi(137ULL),stoi(148ULL)),mkcol2(stoi(26ULL),stoi(145ULL)));
	global_setup.action_4[9] = mkmat2(mkcol2(stoi(95ULL),stoi(238ULL)),mkcol2(stoi(138ULL),stoi(186ULL)));
	global_setup.action_2[10] = mkmat2(mkcol2(stoi(158ULL),stoi(305ULL)),mkcol2(stoi(228ULL),stoi(303ULL)));
	global_setup.action_3[10] = mkmat2(mkcol2(stoi(350ULL),stoi(273ULL)),mkcol2(stoi(380ULL),stoi(112ULL)));
	global_setup.action_4[10] = mkmat2(mkcol2(stoi(169ULL),stoi(307ULL)),mkcol2(stoi(398ULL),stoi(292ULL)));
	global_setup.action_2[11] = mkmat2(mkcol2(stoi(198ULL),stoi(427ULL)),mkcol2(stoi(212ULL),stoi(323ULL)));
	global_setup.action_3[11] = mkmat2(mkcol2(stoi(413ULL),stoi(147ULL)),mkcol2(stoi(207ULL),stoi(109ULL)));
	global_setup.action_4[11] = mkmat2(mkcol2(stoi(317ULL),stoi(104ULL)),mkcol2(stoi(201ULL),stoi(204ULL)));
	global_setup.action_2[12] = mkmat2(mkcol2(stoi(2477ULL),stoi(2023ULL)),mkcol2(stoi(3828ULL),stoi(1446ULL)));
	global_setup.action_3[12] = mkmat2(mkcol2(stoi(1059ULL),stoi(2707ULL)),mkcol2(stoi(502ULL),stoi(2865ULL)));
	global_setup.action_4[12] = mkmat2(mkcol2(stoi(2068ULL),stoi(2456ULL)),mkcol2(stoi(1530ULL),stoi(1855ULL)));
	global_setup.action_twist_2[0] = mkmat2(mkcol2(strtoi("7162594516453352905262709488208"),strtoi("7011190501393501994132146041935")),mkcol2(strtoi("9191917704942915738942923811037"),strtoi("3138456944424184548710837779635")));
	global_setup.action_twist_3[0] = mkmat2(mkcol2(strtoi("4139420820192497085813244422771"),strtoi("7967053796494599236293564809974")),mkcol2(strtoi("4392966125320486618876362200315"),strtoi("6161630640685040368160302845073")));
	global_setup.action_twist_4[0] = mkmat2(mkcol2(strtoi("5448208594922421066181565210788"),strtoi("777225683462809465405647334838")),mkcol2(strtoi("10286881196292221841921057861114"),strtoi("4852842865955116387791982057055")));
	global_setup.action_twist_2[1] = mkmat2(mkcol2(stoi(12ULL),stoi(4ULL)),mkcol2(stoi(6ULL),stoi(1ULL)));
	global_setup.action_twist_3[1] = mkmat2(mkcol2(stoi(4ULL),stoi(5ULL)),mkcol2(stoi(4ULL),stoi(10ULL)));
	global_setup.action_twist_4[1] = mkmat2(mkcol2(stoi(12ULL),stoi(9ULL)),mkcol2(stoi(2ULL),stoi(1ULL)));
	global_setup.action_twist_2[2] = mkmat2(mkcol2(stoi(9ULL),stoi(6ULL)),mkcol2(stoi(9ULL),stoi(8ULL)));
	global_setup.action_twist_3[2] = mkmat2(mkcol2(stoi(1ULL),stoi(14ULL)),mkcol2(stoi(3ULL),stoi(0ULL)));
	global_setup.action_twist_4[2] = mkmat2(mkcol2(stoi(10ULL),stoi(7ULL)),mkcol2(stoi(16ULL),stoi(7ULL)));
	global_setup.action_twist_2[3] = mkmat2(mkcol2(stoi(15ULL),stoi(11ULL)),mkcol2(stoi(42ULL),stoi(28ULL)));
	global_setup.action_twist_3[3] = mkmat2(mkcol2(stoi(18ULL),stoi(15ULL)),mkcol2(stoi(24ULL),stoi(26ULL)));
	global_setup.action_twist_4[3] = mkmat2(mkcol2(stoi(18ULL),stoi(38ULL)),mkcol2(stoi(9ULL),stoi(25ULL)));
	global_setup.action_twist_2[4] = mkmat2(mkcol2(stoi(28ULL),stoi(26ULL)),mkcol2(stoi(64ULL),stoi(51ULL)));
	global_setup.action_twist_3[4] = mkmat2(mkcol2(stoi(48ULL),stoi(4ULL)),mkcol2(stoi(58ULL),stoi(32ULL)));
	global_setup.action_twist_4[4] = mkmat2(mkcol2(stoi(8ULL),stoi(75ULL)),mkcol2(stoi(26ULL),stoi(71ULL)));
	global_setup.action_twist_2[5] = mkmat2(mkcol2(stoi(4ULL),stoi(89ULL)),mkcol2(stoi(118ULL),stoi(153ULL)));
	global_setup.action_twist_3[5] = mkmat2(mkcol2(stoi(117ULL),stoi(28ULL)),mkcol2(stoi(90ULL),stoi(41ULL)));
	global_setup.action_twist_4[5] = mkmat2(mkcol2(stoi(0ULL),stoi(150ULL)),mkcol2(stoi(101ULL),stoi(0ULL)));
	global_setup.action_twist_2[6] = mkmat2(mkcol2(stoi(230ULL),stoi(68ULL)),mkcol2(stoi(48ULL),stoi(9ULL)));
	global_setup.action_twist_3[6] = mkmat2(mkcol2(stoi(169ULL),stoi(61ULL)),mkcol2(stoi(230ULL),stoi(71ULL)));
	global_setup.action_twist_4[6] = mkmat2(mkcol2(stoi(18ULL),stoi(216ULL)),mkcol2(stoi(144ULL),stoi(221ULL)));
	global_setup.action_twist_2[7] = mkmat2(mkcol2(stoi(121ULL),stoi(162ULL)),mkcol2(stoi(194ULL),stoi(150ULL)));
	global_setup.action_twist_3[7] = mkmat2(mkcol2(stoi(137ULL),stoi(217ULL)),mkcol2(stoi(197ULL),stoi(135ULL)));
	global_setup.action_twist_4[7] = mkmat2(mkcol2(stoi(253ULL),stoi(160ULL)),mkcol2(stoi(31ULL),stoi(18ULL)));
	global_setup.action_twist_2[8] = mkmat2(mkcol2(stoi(123ULL),stoi(180ULL)),mkcol2(stoi(26ULL),stoi(160ULL)));
	global_setup.action_twist_3[8] = mkmat2(mkcol2(stoi(221ULL),stoi(140ULL)),mkcol2(stoi(246ULL),stoi(63ULL)));
	global_setup.action_twist_4[8] = mkmat2(mkcol2(stoi(147ULL),stoi(260ULL)),mkcol2(stoi(109ULL),stoi(136ULL)));
	global_setup.action_twist_2[9] = mkmat2(mkcol2(stoi(102ULL),stoi(165ULL)),mkcol2(stoi(123ULL),stoi(205ULL)));
	global_setup.action_twist_3[9] = mkmat2(mkcol2(stoi(15ULL),stoi(293ULL)),mkcol2(stoi(26ULL),stoi(293ULL)));
	global_setup.action_twist_4[9] = mkmat2(mkcol2(stoi(294ULL),stoi(253ULL)),mkcol2(stoi(114ULL),stoi(13ULL)));
	global_setup.action_twist_2[10] = mkmat2(mkcol2(stoi(134ULL),stoi(14ULL)),mkcol2(stoi(527ULL),stoi(429ULL)));
	global_setup.action_twist_3[10] = mkmat2(mkcol2(stoi(222ULL),stoi(203ULL)),mkcol2(stoi(62ULL),stoi(342ULL)));
	global_setup.action_twist_4[10] = mkmat2(mkcol2(stoi(214ULL),stoi(462ULL)),mkcol2(stoi(27ULL),stoi(349ULL)));
	global_setup.action_twist_2[11] = mkmat2(mkcol2(stoi(228ULL),stoi(441ULL)),mkcol2(stoi(166ULL),stoi(371ULL)));
	global_setup.action_twist_3[11] = mkmat2(mkcol2(stoi(46ULL),stoi(270ULL)),mkcol2(stoi(324ULL),stoi(554ULL)));
	global_setup.action_twist_4[11] = mkmat2(mkcol2(stoi(28ULL),stoi(384ULL)),mkcol2(stoi(253ULL),stoi(571ULL)));
	global_setup.action_twist_2[12] = mkmat2(mkcol2(stoi(415ULL),stoi(76ULL)),mkcol2(stoi(82ULL),stoi(192ULL)));
	global_setup.action_twist_3[12] = mkmat2(mkcol2(stoi(18ULL),stoi(519ULL)),mkcol2(stoi(407ULL),stoi(590ULL)));
	global_setup.action_twist_4[12] = mkmat2(mkcol2(stoi(161ULL),stoi(429ULL)),mkcol2(stoi(103ULL),stoi(446ULL)));
	global_setup.action_twist_2[13] = mkmat2(mkcol2(stoi(59ULL),stoi(187ULL)),mkcol2(stoi(607ULL),stoi(560ULL)));
	global_setup.action_twist_3[13] = mkmat2(mkcol2(stoi(152ULL),stoi(69ULL)),mkcol2(stoi(551ULL),stoi(468ULL)));
	global_setup.action_twist_4[13] = mkmat2(mkcol2(stoi(585ULL),stoi(594ULL)),mkcol2(stoi(331ULL),stoi(34ULL)));
	global_setup.action_twist_2[14] = mkmat2(mkcol2(stoi(229ULL),stoi(385ULL)),mkcol2(stoi(711ULL),stoi(514ULL)));
	global_setup.action_twist_3[14] = mkmat2(mkcol2(stoi(209ULL),stoi(719ULL)),mkcol2(stoi(217ULL),stoi(535ULL)));
	global_setup.action_twist_4[14] = mkmat2(mkcol2(stoi(638ULL),stoi(612ULL)),mkcol2(stoi(87ULL),stoi(105ULL)));
	global_setup.action_twist_2[15] = mkmat2(mkcol2(stoi(721ULL),stoi(620ULL)),mkcol2(stoi(290ULL),stoi(106ULL)));
	global_setup.action_twist_3[15] = mkmat2(mkcol2(stoi(595ULL),stoi(756ULL)),mkcol2(stoi(674ULL),stoi(233ULL)));
	global_setup.action_twist_4[15] = mkmat2(mkcol2(stoi(27ULL),stoi(645ULL)),mkcol2(stoi(29ULL),stoi(800ULL)));
	global_setup.action_twist_2[16] = mkmat2(mkcol2(stoi(102ULL),stoi(413ULL)),mkcol2(stoi(631ULL),stoi(839ULL)));
	global_setup.action_twist_3[16] = mkmat2(mkcol2(stoi(523ULL),stoi(442ULL)),mkcol2(stoi(780ULL),stoi(419ULL)));
	global_setup.action_twist_4[16] = mkmat2(mkcol2(stoi(27ULL),stoi(760ULL)),mkcol2(stoi(147ULL),stoi(914ULL)));
	global_setup.action_twist_2[17] = mkmat2(mkcol2(stoi(1194ULL),stoi(2293ULL)),mkcol2(stoi(1615ULL),stoi(1163ULL)));
	global_setup.action_twist_3[17] = mkmat2(mkcol2(stoi(1893ULL),stoi(2199ULL)),mkcol2(stoi(2104ULL),stoi(465ULL)));
	global_setup.action_twist_4[17] = mkmat2(mkcol2(stoi(1929ULL),stoi(789ULL)),mkcol2(stoi(552ULL),stoi(428ULL)));
	global_setup.action_two_2 = mkmat2(mkcol2(strtoi("28623943409757479661"),strtoi("16627037906587991710")),mkcol2(strtoi("30730418368513759301"),stoi(8269544737661623571ULL)));
	global_setup.action_two_3 = mkmat2(mkcol2(strtoi("28623943409757479661"),strtoi("16627037906587991710")),mkcol2(strtoi("30730418368513759301"),stoi(8269544737661623571ULL)));
	global_setup.action_two_4 = mkmat2(mkcol2(strtoi("28623943409757479661"),strtoi("16627037906587991710")),mkcol2(strtoi("30730418368513759301"),stoi(8269544737661623571ULL)));
	for (int i = 0; i < 3; ++i) {
		fp_enc( &(&(&torsion_basis_two[i])->x)->re, &(torsion_basis_two_uintbig[i][0][0]) );
		fp_enc( &(&(&torsion_basis_two[i])->x)->im, &(torsion_basis_two_uintbig[i][0][1]) );
		fp_enc( &(&(&torsion_basis_two[i])->z)->re, &(torsion_basis_two_uintbig[i][1][0]) );
		fp_enc( &(&(&torsion_basis_two[i])->z)->im, &(torsion_basis_two_uintbig[i][1][1]) );
	}
	for (int i = 0; i < 3; ++i) {
		fp_enc( &(&(&torsion_basis_sum[i])->x)->re, &(torsion_basis_sum_uintbig[i][0][0]) );
		fp_enc( &(&(&torsion_basis_sum[i])->x)->im, &(torsion_basis_sum_uintbig[i][0][1]) );
		fp_enc( &(&(&torsion_basis_sum[i])->z)->re, &(torsion_basis_sum_uintbig[i][1][0]) );
		fp_enc( &(&(&torsion_basis_sum[i])->z)->im, &(torsion_basis_sum_uintbig[i][1][1]) );
	}
	for (int i = 0; i < 3; ++i) {
		fp_enc( &(&(&torsion_basis_twist_sum[i])->x)->re, &(torsion_basis_twist_sum_uintbig[i][0][0]) );
		fp_enc( &(&(&torsion_basis_twist_sum[i])->x)->im, &(torsion_basis_twist_sum_uintbig[i][0][1]) );
		fp_enc( &(&(&torsion_basis_twist_sum[i])->z)->re, &(torsion_basis_twist_sum_uintbig[i][1][0]) );
		fp_enc( &(&(&torsion_basis_twist_sum[i])->z)->im, &(torsion_basis_twist_sum_uintbig[i][1][1]) );
	}
	for (int j=0;j<13;j++){
		for (int i = 0; i < 3; ++i) {
			fp_enc( &(&(&torsion_basis[j][i])->x)->re, &(torsion_basis_uintbig[j][i][0][0]) );
			fp_enc( &(&(&torsion_basis[j][i])->x)->im, &(torsion_basis_uintbig[j][i][0][1]) );
			fp_enc( &(&(&torsion_basis[j][i])->z)->re, &(torsion_basis_uintbig[j][i][1][0]) );
			fp_enc( &(&(&torsion_basis[j][i])->z)->im, &(torsion_basis_uintbig[j][i][1][1]) );
		}
	}
	for (int j=0;j<18;j++){
		for (int i = 0; i < 3; ++i) {
			fp_enc( &(&(&torsion_basis_twist[j][i])->x)->re, &(torsion_basis_twist_uintbig[j][i][0][0]) );
			fp_enc( &(&(&torsion_basis_twist[j][i])->x)->im, &(torsion_basis_twist_uintbig[j][i][0][1]) );
			fp_enc( &(&(&torsion_basis_twist[j][i])->z)->re, &(torsion_basis_twist_uintbig[j][i][1][0]) );
			fp_enc( &(&(&torsion_basis_twist[j][i])->z)->im, &(torsion_basis_twist_uintbig[j][i][1][1]) );
		}
	}
    	for (int i=0;i<3;i++){
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->x)->re, &(torsion_basis_ted_sum_uintbig[i][0][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->x)->im, &(torsion_basis_ted_sum_uintbig[i][0][1]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->y)->re, &(torsion_basis_ted_sum_uintbig[i][1][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->y)->im, &(torsion_basis_ted_sum_uintbig[i][1][1]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->z)->re, &(torsion_basis_ted_sum_uintbig[i][2][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->z)->im, &(torsion_basis_ted_sum_uintbig[i][2][1]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->t)->re, &(torsion_basis_ted_sum_uintbig[i][3][0]) );
    		fp_enc( &(&(&torsion_basis_ted_sum[i])->t)->im, &(torsion_basis_ted_sum_uintbig[i][3][1]) );
  	}
  	for (int i=0;i<3;i++){
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->re, &(torsion_basis_twist_ted_sum_uintbig[i][0][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->x)->im, &(torsion_basis_twist_ted_sum_uintbig[i][0][1]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->re, &(torsion_basis_twist_ted_sum_uintbig[i][1][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->y)->im, &(torsion_basis_twist_ted_sum_uintbig[i][1][1]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->re, &(torsion_basis_twist_ted_sum_uintbig[i][2][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->z)->im, &(torsion_basis_twist_ted_sum_uintbig[i][2][1]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->re, &(torsion_basis_twist_ted_sum_uintbig[i][3][0]) );
    		fp_enc( &(&(&torsion_basis_twist_ted_sum[i])->t)->im, &(torsion_basis_twist_ted_sum_uintbig[i][3][1]) );
  	}

}

