// p = 73743043621499797449074820543863456997944695372324032511999999999999999999999
//
// p-1 = 2 * 3^53 * 43 * 103^2 * 109 * 199 * 227 * 419 * 491 * 569 * 631 * 677 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 4283
// p+1 = 2^33 * 5^21 * 7^2 * 11 * 31 * 83 * 107 * 137 * 751 * 827 * 3691 * 4019 * 6983 * 517434778561 * 26602537156291

#include "constants.h"

const long class_mod_4 = 3;
const long two_tors_height = 33;
const bool need_even_commit = false;
const bool need_even_chall = false;
const long security_level = 128;

const long signing_length=1000 ;

const long signing_length_two_tors_height_step = 31;
const long last_step_length = 10;

//should be two_tors_height + delta
//but is a bit tricky right now
const long len_tail = 10;

const char* p_str =
  "73743043621499797449074820543863456997944695372324032511999999999999999999999";
const char* all_the_torsion_str =
  "5438036482562421961818827832407841404369503851311267868831048027176534732013273363338731215746134716355745303230004110609255351934976000000000000000000000";

const uintbig p_plus_odd_cofactor = { 0x68cd740600000000, 0x0016c5bcbd22f015, 0, 0 };
const uintbig p_minus_odd_cofactor = { 2, 0, 0, 0 };
const uintbig p_even_cofactor = { 0xa52ca964a8652149, 0x1bb9479de8d8027c,
				  0xdb3c54c8592e3b52, 0x51848ab2 };
const uintbig p_even_cofactor_minus_1 = { 0xa52ca964a8652148, 0x1bb9479de8d8027c,
				  0xdb3c54c8592e3b52, 0x51848ab2 };
const int p_even_cofactor_minus_1_limb = 3;
const uint64_t p_even_cofactor_minus_1_highest = 0x40000000;

const uintbig Pxr = {0x61d4f73ae4c50a7a, 0x94afdf85986ab52f, 0x81ddac60a790f21a, 0x88555095527404f7};
const uintbig Pxi = {0,0,0,0};
const uintbig Pyr = {0,0,0,0}; 
const uintbig Pyi = {0xc3e0d9206ba0f388, 0xd1942b76fe67a141, 0xd6a7b5fe81127233, 0x6227a343e295edd9};

#define M_LEN 19
const long p_minus_len = M_LEN;
const long p_minus_fact[M_LEN] =
  {  3, 43, 103, 109, 199, 227, 419, 491, 569, 631, 677, 857, 859, 883,
    1019, 1171, 1879, 2713, 4283 };
const long p_minus_mult[M_LEN] =
  { 53,  1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
       1,    1,    1,    1,    1 };




#define P_LEN 12
const long p_plus_len = P_LEN;
const long p_plus_fact[P_LEN] =
  {  5, 7, 11, 31, 83, 107, 137, 751, 827, 3691, 4019, 6983 };
const long p_plus_mult[P_LEN] =
  { 21, 2,  1,  1,  1,   1,   1,   1,   1,    1,    1,    1 };


// the multiplicities to take to obtain log2(p) bits of torsion (for commitment)
const long p_minus_mult_com[M_LEN] =
  { 0,  1,   2,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,   1,
       1,    1,    1,    1,    1 };
const long p_plus_mult_com[P_LEN] =
  { 0, 2,  1,  1,  1,   1,   1,   1,   1,    1,    1,    1 };


// the multiplicities to take to obtain log2(p)/2 bits of torsion (for challenge)
const long p_minus_mult_cha[M_LEN] =
  { 53,  0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
       0,    0,    0,    0,    0 };
const long p_plus_mult_cha[P_LEN] =
  { 21, 0,  0,  0,  0,   0,   0,   0,   0,    0,    0,    0 };

const int win = 5;
const int Strategy[5] = {4, 1, 2, 1, 1};

const uintbig Lookupbasere = {0x603e2a275707553c,0x5d1c1b934dd7e707,0x3dd0d46d31bdc9c,0x8f47dc9c611528a5};
const uintbig Lookupbaseim = {0xfb5a0c95da7f9c94,0xe96ea7df88001688,0x5e16ebb480ddbb88,0x83b699cd66ec0f44};
const uintbig Lookuptable[160] ={
{0xe1f344d7acf51499,0x4f30b30ea8a1fe1a,0x4a03f03b8ded72b1,0x1043b7ee71bc8022},
{0xacf511077889bd1f,0xed5a9b3f51e9f0bc,0x2af4a239b9ba457d,0x2ca9029407de1d4d},
{0x2a174d0f2e226f3d,0xff8b3b8b73eb34af,0xcb9b8b9b0b928fc0,0x2a31fbb8c2cbfbc4},
{0xb1ff445c528259fb,0x308148672dcfb69a,0x570169dd358ac350,0x63a1d52e8bbcf0e6},
{0xb5e718767921f7d1,0x41acc9780f7a1808,0xdf81d71db9f2d3c7,0x754ed2c3e5ab6719},
{0xd82dc0789e9c53cd,0xc9a697e752be1ed6,0x9f5976ec1f0bbd05,0x580f521885ec1f1f},
{0x9fd7e32fc86c7f55,0x32ea24408bba6f7,0x18920509c9482d88,0x7bea2548e83a628f},
{0xada2c50c2e3f6e61,0xa4010586a94ed01b,0xda05fcfdc724b156,0x4a074676b1858404},
{0x8c294ee700066e8,0x4caaea5001cc999a,0xb67a730a713f5beb,0x5e1e0c8e1f166846},
{0xc5006b2c2034a685,0x9b5a114ff94c8ffa,0xe44531f9cf4ebd16,0x157ebb16cbb2c98e},
{0xa819595f0056431a,0x6caf1202666aa152,0xe2aafa51ca62c2b8,0x2d1bac17a5dfa0d2},
{0x6a32482d21e418ab,0xd664782a72a51612,0xb6dc590a85f1768a,0x50eeae43c85c91a0},
{0x4c5112fa37360d69,0xeb5bb699d8a07a9d,0xb5e7136ae2158213,0x7bfd0a902f8f633e},
{0xa34426b34fbd134c,0xeb85acc67e2a79c5,0x25f34666fa01fa38,0x52c8097e35f61bd9},
{0x4c0d8f8b9ec98ac6,0xc62afed889f17778,0xf027ec187a49f1be,0x4e19b4bac7f8927b},
{0xfda38737c4c90a6,0x9631a39d7288ecaa,0x6b9698a4190ef8f4,0x1ccd338f86074ed8},
{0x97c1cf9f5cec2d0,0x2521ff5a6cc7995c,0x61abd4120b07605a,0x62319d6185a9116a},
{0x62041dd215bf5e3e,0xbd9c46965260aef5,0x59d60c6a3f8e2994,0x17a77b3e2d657e17},
{0x6f5d104991136fa8,0xfcd9ca7e42e28e94,0x2cf09d808c2c031b,0x7100e6103d1b96ed},
{0xa919471d6c95ec1b,0xc34744f2dab83080,0xa1bae5dd95553b34,0x14b8d4e6cdf35bdd},
{0xfdf1ef7e9207d6ca,0x13ba561ed50592a9,0x5398dd43f5277dce,0xa0ea574d06dd59f5},
{0x767fd90fba8afa9c,0x6c9c9c367b71b68f,0xeb30c0b7693b14f6,0x12da805f6260821e},
{0x55703185b2c38c04,0x5cad26da109bb44d,0x57e4bbc494ce07c6,0x5c85e41aa74f8b2b},
{0x7bc0cc68f2e2620a,0x689eb93ececbef69,0xb0ed59970df281b6,0xd09e40c246aa9d8},
{0xe7ce6fc98091a622,0x473c5b1f9acfe93f,0xc271c1daba847ae6,0x85fa06cd11329b3c},
{0x958e21eca9ccfdf0,0x7efcee83025954f1,0x5fe4b1776811927,0xf7b103d957a08a8},
{0xc42f55a9077eebb6,0x6acf3cbb093f6f90,0x90fe509d74820801,0x9451b87a78505966},
{0x9307577b14ddcd4c,0xb752811cb29f27d2,0x685382f3f44b3f04,0x7cd869a311798c57},
{0x48fc294c1ca87372,0x5cc31943fd87bd2a,0xa4829c1e0dc9e573,0x506ae65e5c801516},
{0xcf04012e0acf4dd2,0x490362d13a8e8372,0xc17b2f1c2e952032,0x9b0dc97e35bcfb01},
{0x72100f5aa46158ac,0xbbd5279328a675d8,0xf388b6b63d0d49e2,0x879b9351d32fc82a},
{0xae4c9c9c3e0bb70d,0x8f6c52cfb1bd8552,0xf8683edb74827be3,0x14357d03601104f6},
{0xa85df2c7bbbae7ae,0xcd4dac0db5efe93a,0x6be80266f9ea699e,0x5df535940fc2c1d6},
{0x1b747afd7aded50,0x6081f06c953df881,0xbd53e8faac87f620,0x9ece5f7748973795},
{0x8440d9a262b35fe2,0x3c4fd164e2987d16,0xfb9fecc29a9f75aa,0x77f80a0a963ccf4c},
{0xe317f8a1d3b8efec,0x98a135f2954e295d,0x805bac073a31571f,0x28fda83f06a82178},
{0x4390a6826b02123,0x6615c9a3a5f9ad8e,0xbd4d00bfe6453dd,0x60751a68aacc19db},
{0x26e309cb240ed84,0x927717cc23bd0595,0x3d2735b935c557d3,0x99e41036823b6d7f},
{0xe85c6789da4b69f7,0x62b3b712343367d6,0xdf189e9d4becb023,0x72c4f4ac68761762},
{0x1acd53da2a037080,0x54afbb4900880762,0x8dd76593ca26485c,0x4e08de1ec212cfe6},
{0x6438f3d67718a4ca,0xfa528c9e0dc37bec,0x7216856221298dc9,0x50be72cbc45df871},
{0x1b334fb3fd43deee,0xa41472a19b8b5dd8,0xe0fc60f3efa15f54,0x575a3150e0a16504},
{0x3f78349e3088c377,0x9b54d48252009390,0xb2d8d00330a0fa3,0x9ea505bfb464daad},
{0xc3dfaf4bcbf9a4e2,0xeedaa4c553c9cf40,0x1b17d0cd3080bb0,0x1bfd702586120d54},
{0x27cd65b5e49a37c7,0x795b4a21e27c3827,0xabe080f79a0605a0,0x5fbe4d07121f8a30},
{0xb2166fa38a06da87,0xa2eeec1048c1ff60,0xdff27d902f42c4af,0x73192c038d0a6d3d},
{0x53543b54176200e4,0x3e71166a66eb856a,0xab0c2a2269b047c0,0x5b6e4fbf7b797564},
{0x48dbeb9e86290634,0xa66207c36d4a2f83,0xba4da6df3e8e2ef,0x47d29ac89137e0cf},
{0x1a8aaeb0ae5770a6,0xffa043ca37828501,0x6b86d94c0be06757,0x1a58a26466e3f9ff},
{0xcaa78d927e2c1f93,0x6a2e78665aa9b416,0x702419f137bdc9bd,0x8c4ec864c2ec96fb},
{0x7c858be2db29bc93,0x6400f412e1c8fb88,0xf68af239b4c7851e,0x80ca53a2cdf3b1d5},
{0x9c2ec7d023e6e2ec,0x654dc245bd2876a0,0x9e7e5b24ccc89680,0x82e1699379cecf56},
{0xc857a80084a3236d,0xcfe508e4f369f3e3,0xcaca0089d1a544d1,0x9ea0f11c5b511a0c},
{0x8966b1b8dc606564,0x46678d6546df2c39,0x3cc1fb5b5470826,0x757a3227f202b716},
{0x70cc48ce8c1b5709,0xb49becf58449195,0x6f282af128ba9586,0x6c728dbc9e29f4c7},
{0xbb92b55c1b116a95,0xd08c280a0f608869,0x1ed3740deb1bbe0b,0x2533dd56e6646b3f},
{0xee0170d7fcc389a9,0xaaebf4bef1d4fd86,0x577ca00361ca672,0x6d3d25963be4734e},
{0x45a60835953148e7,0xf76c330eb441e4d7,0x3abda0b5546cc5f8,0x883336a21785d7},
{0x51d319d797c8bf5b,0x89e3debde0e27c3b,0x7b8d1fb274c5dd5f,0x5843793833dc6b33},
{0xf4b2fc64ce11f9e,0x11e41a7e9744f53c,0xcabee9f5e6a44074,0x501171eb96b6de93},
{0x8e4d106d5dd9d795,0x7e50b6b3388106fa,0xe081957b129afe88,0x9be204d01e18c46a},
{0xfa06b2a49f99d9c,0xf6416101fd9ad79c,0xb99bbe66610e7508,0x34ea92f9fbd9f2fe},
{0x39205f26132e26dd,0x349212bf69adb4a7,0x26e82b02025bc6d8,0x5d271689cd95862e},
{0x1155eef06effa820,0xff16967ec4cc2f22,0x161a3e72bdb6a6ad,0x2d7bb9900566e83b},
{0xadc418fdfc448c48,0x9029c23d8a267d42,0x85fec463209bbfe6,0x91aa78e09d9828e6},
{0xbd25ef5816150a7f,0x2130738e097aecf4,0xe1238b399fcd2687,0x15d91b62a63cb204},
{0x5faea1cba7714a62,0xfffc2cfaaeb8d556,0xdd03edd188a872ee,0x23ed688d65aaa48a},
{0xfd64980657e16a0,0x991be20faef485de,0xee8f428ca81da209,0x25701e5df4a06cb8},
{0x58189029ee7d2cba,0x35aff62a78f8b150,0x9c21d38cc21ac101,0x3bea57ad0c477283},
{0x1d16e92340401fa1,0xffda9a706570b667,0xbd6d41fde37ea51a,0x6a21430c964c1994},
{0xaecaec64514d920f,0x9f88d58f2272e00,0x9cf83152315a6cd0,0x76f1fdc122557ae8},
{0x3134139468115940,0xcc99a88bdc6d77e2,0x8053f8206aca54ac,0x392bb28d53575dea},
{0x4d50a2c45d66324d,0xab7655997d919b7c,0x356031c7b27dabf1,0x97073b524ca89928},
{0x97b4e8ef97159df7,0xd9ebeb20af041359,0x6151fb4617791884,0x9e37b9de45d09474},
{0x172bb8667b7eb564,0xc29e015cb26e12ea,0x8f513671026a8dea,0x32007a1792545680},
{0x6e5766ee4a5109a2,0x6e98f566a4aa6a7c,0xf0b40f709a759d55,0x89d634599272c223},
{0x7aba610d5f69cf5a,0x8864c80db602f6f2,0x79c94599b102d64f,0x3fe4a849a55603c4},
{0xe81579e324a14887,0x1b2838def9c09fbe,0x99bbf6b18f623e06,0x70effa100a0070f7},
{0xed3e8080bad3000d,0xcc7d167eef23d98b,0x477a92f9e9a12391,0x99f354b6357760f1},
{0x3fced90b811ce9eb,0x2b7b81a765cbdec6,0xa7b9e51cb4cd45a7,0xa6decc6f00e623a},
{0x10599183354bd3c0,0xff10158b93354988,0x7a703772d81ce97,0x3d5b33ff3f987d4e},
{0xc72fd2e47d529a71,0x8e272800c856e0d8,0x8d5d17cef267b8ea,0x55e9256d9bb2de84},
{0xdb344f9997326a36,0x328118b614d7e521,0xe699be87a523b21c,0x9f65f81fd5a0bb66},
{0x64b7623f693b0d64,0x9365cc3fbeec5caa,0x903792ac4ef7dee8,0x202d366a834aad83},
{0xecb19e7ab97a0555,0xb36551ab1fd9563c,0xee18f550867fea2e,0x4169112f3b9fb674},
{0xa4f587b324241f5a,0x7a912803c02f4afa,0xd61fe5822681c95f,0xcc490394b9e15a3},
{0x325cff9a1440c1e6,0x134caf8409f5d47a,0x732f5c3378d77eab,0x71587796d3cf0e3c},
{0xfbc9c42ce178acd5,0x759f72497ecd71cb,0x3ab490b13fb9390c,0x9c1060fe5cd917de},
{0xac6f7ff5c75a4d48,0x27bf6f633940e107,0x1299e4bc868d9f44,0x83351d827ef76935},
{0xf74968089bfe2f54,0xfe6ad0c1463c5911,0xb34d0d813fd055a7,0x5f675fb2f0e9d74},
{0x1d719a5beb09d9dc,0x7eab085c94d46631,0x1d0cb6ca9c8f7938,0x93a06d135c8cce1e},
{0x1ecffc11ad587c2e,0x5f8ade4f86cf9732,0xf3c422dc01a33121,0x135281c6fce49b31},
{0x3a015b3cb8f9cb7d,0xe8f3002bb36ad663,0x63f75c555e9086a8,0x46e5acbdf22563dd},
{0x9097adc5b33a08db,0xa5ce9022a1684352,0xca494496c71253d,0x9a0b8b0c82b3ca54},
{0x690065de6cd792dd,0x6249537c2a4cf3ef,0x93260928262e9a39,0x40a408c971192e01},
{0x822f0ce93d3c9e5b,0x1151a70efdab3c74,0xad1d923efb51b8de,0x1fd85e16fcf67265},
{0xa32a27b05dfa645f,0xf641610f9a82852d,0x4c2f491808dfb3f,0x7b9f90b13e3f55db},
{0xcb3fe27830aec698,0x8c4612f677762858,0x7770f0d3d0fc045b,0x1fb58d04ea54576a},
{0x6d72e0f4b3c5ec04,0xf1118c3ca8ffdb0b,0xd955e0823aee083f,0x5c4326105e962ade},
{0x2ca88eec24213501,0x45ec676201b75484,0x7e378f4ed06ed8b1,0x7b76034bd2e25d3},
{0x6273aaf690159dc2,0x248ba7bdf8f78ba1,0xaf4f9cb4ed77d031,0x26ca6166a2c58347},
{0xb321c764d1647e23,0x7880ecb38dffb7bd,0x55498ebf38c0b807,0x5c2e2aacbb133cde},
{0xe71de363fb339b44,0x4f618480eefbac96,0xcde93cf9ec8d45e6,0xa19fac0cc681c1b9},
{0x225f6a9b3ad7d0c7,0xaa36e33f612c6f3b,0xba53d0ab68d42e49,0x1d971928814498a5},
{0x4749c24ec4797918,0x2381e247f14644b8,0x8453d587919e531,0x7aacafe376d1f7a5},
{0x58aea7e313387978,0x979e302fac354e3c,0xebad6da5ab82d45a,0x347d26fb7832901e},
{0x41a6b7be1690cd5d,0xbd1923f379cbfcbc,0x3846ad68519a3770,0x2417d44ef7072d9a},
{0x609a7ae58e2b9639,0x7e7a576e6a1c6f4c,0x9d804bf3634a1821,0x4bcc2448c9c556b0},
{0xf0607a3425d004ac,0xde7c2a20f13a6f19,0x714afc41aa4b1607,0x8abe78ee154d7704},
{0x5e01fc74261f1c57,0xc15cd148ad766a48,0xbb91e5ee3671cc36,0x26cdbaa30336db80},
{0xf48bb0e334eb5cb8,0x3995a9cb9324afaa,0xd8024a221fca8b94,0x38aa2fd28c187fd2},
{0x633b0e0709748aa0,0x22f07180339f9ed2,0x9e1e17a3a901aef,0x18159376160e2b7f},
{0x4f48bfc92bef1c38,0x778166c9378ddec3,0xa504fd9904fa81c4,0x11c73f3d89b0a318},
{0x82b787faac5c41d6,0x435eb68d149793a9,0xeaf4ec148a063ab8,0x79eb5186d6d58c86},
{0x523b14983b976e37,0x1c4e1842cf45b6ce,0xf7696e373f69a6ba,0x19b551361b06e39b},
{0x646868ece0d937a7,0xd5a8f9c46e9b3260,0x3b0f07058c8e434e,0x38498513b9087724},
{0xcff3d359fa2abe7a,0xd3fe2b0e158afde0,0x37b7d0652ff3d191,0x3e4ae0c736fc64f2},
{0xf6f0ff2b841bd663,0x6e55e22525a39fe5,0x3afb1744ea3e295a,0xa265768064b1b8c9},
{0x4a8111ab5c0f0e35,0x7df25e41d1f4b781,0x66a61e67a69a14fd,0x8a0862eb7d31ae61},
{0x3f2a91a5234f7dfd,0xc7d18d71893434f4,0x4db5b40d215b7528,0x81d78f6ccd84b465},
{0x584ea6f63313ca1,0x184be1a4456ca812,0x96efd01e8c6ad529,0x3f57ad3a07ddc960},
{0x76bdc172def7ff18,0x15909a8ef43ca33e,0x9348fca86076b53f,0x6d9723f8fd4abe00},
{0x588bc4faf8f1c358,0xfcc2fa01084a2d2,0x8cdc203e38066ab7,0x925fcb0ac4e54bd2},
{0xe0eccef03200119b,0x560735dfecb98102,0xebf6d5ea36d37b87,0x1aff8a1baecd40ae},
{0x3b802b0d7cea41f2,0x480864e92fc23178,0xee6a35b5105c8959,0x316597cde9f760b8},
{0xd3164d05bcc1a1e6,0x45963bfb3fb5fedf,0x8905daee922beec3,0xa422913909a551c},
{0xac09bb26c8d71768,0x4d24d0aa9415f6a6,0xed8347e3293063e2,0x850e4185d6f82f38},
{0x21dbfe3f3351e430,0xd6151178ddb7e645,0xd4402df81109e6e9,0x3f881692e983120f},
{0x821d4d7560470ce1,0x8166c4b245121c1e,0x5e91edf8b38d0da1,0x7cc437debd61a1b},
{0x688887cb522a0b5f,0xead88c76dd1fe9d6,0xf7415f64113212cd,0x101e9dc329f4f487},
{0xf9c9e56a079052f9,0x31776a33095f32b0,0xf11837c4bb1d5ba0,0x69e4e72643b3dab5},
{0xb4a83ae8b192e3cd,0x1a016164acadf1eb,0x5e468b89d941011d,0x86582bc20ca88bcb},
{0xd83e28c9c18e0c78,0x7e7931812393ae52,0x4522fb5a7a2ffcda,0x766df7b832c84f87},
{0x92427208cfb1d8f8,0xdc3d73c6f7044eea,0x9c615ab5e82fbb67,0x21b8b267ad4205b},
{0x3139d913576c8633,0x1c8e36b9bd61c19b,0xa4e882e260fe76d6,0x782d640602e8e76},
{0xea6a654f25d6f8f9,0xcd311ecce22abfbc,0x9c686148f24ebaa2,0x636c5038781574bd},
{0xba3f842c03775d92,0xd6f7d53650fe4333,0xfb47948e926b7b73,0x64054c0420ec54a6},
{0xcc53b18c59209237,0xca62c21157979b69,0x7ee6910ae9093a56,0x29922304285a3434},
{0x67c4523795a8cd20,0x7a122c4dc005d529,0x8392f349232c4122,0xb742ec72a258664},
{0x42bec21eb59fa192,0x64a74f13d039d6bb,0x54f9d452bac11c45,0x373e0f03f54c0613},
{0x48f9f4a11a964a71,0x2b9c1aa68c96c6d1,0xdff48f1427d68cac,0x3f4322ad1a743995},
{0xa4db70f2e78b1729,0x84da9696e3013935,0x3153d83b29a6ce57,0x67e4585341f3cf6},
{0xea8dba3009e162f,0x69c325a99a99a373,0x12e6388189c391cf,0x4fb7326980c7bd0a},
{0xef7d7ce41f263ea8,0xb2f863bd112efb02,0x872ecd47c30490b,0x943746d0416888c9},
{0xa0e898a3d16b4f63,0x1e99611bd49bd166,0x6f95f9f20decd637,0x80a911ba6bd490e5},
{0xa3e224d4b0f1ee67,0xc48806c206858dd5,0xcd492307d148e182,0x4b8870ce01d951b3},
{0x52584089718edcb3,0xc8eaac761c364ce2,0x9b309c37304179f2,0x243d73acc6672bc3},
{0x37a373323b214533,0x583f07a1e9834f23,0x7e45a58b02fe54bc,0x1ac50f9fe33b476e},
{0xd1c50051685e0b70,0xa539d15da4255e23,0x57687764907a9f8d,0x5c29bdc45592225f},
{0x5bef8541b37f0b5d,0x53faafb9eecbf8c2,0xf121b4d8da5eadc4,0xa19e6ce660a1a085},
{0xe5e99de49f93682c,0xfd84bb6169ddb4e8,0x641c1abf3a906e9f,0x407192083098c36c},
{0xefc62ead91018752,0xddcc9b74d50fa313,0x4bfc539ac290213d,0x76fe2be060df86fb},
{0xeb940d15d8cf5a49,0xaffe00b4547ed0dd,0x8ddf9992925b601c,0xe7baac484d9261b},
{0x9c7d289089ce3721,0x8e631c234a8c68e0,0x4faa925aa714b557,0x38bb87a1a7cfc51b},
{0xa7b3c63de529d852,0xb37889546fca511d,0x72355dbae8f84753,0x955d33bacbe86eba},
{0x5f6c4d65fd1deb40,0xfc4a7d0ec77e53e8,0x3782a7a67a793442,0x80b260600d9e3aa},
{0x4736ea419e98366d,0xda4f2bea94f90ef2,0xcc056d2c28817bd3,0x97d359dc16380d1},
{0xb9a23aeaaf6387fc,0x98960c1dab033fea,0x9ed08e0fec85afba,0x13782228752adf40},
{0xf853055cbb5640d4,0x7c85956c58bb67c2,0xb88a1ccfb90e8849,0x205f8d768f6832ad},
{0x7db0022b4df904d4,0x90f64670d73d3780,0xdb7e4cef28a58d5c,0x9f81e2f8c72b7ffb}
};
const uintbig LookuptableLR[32] = 
{
{0xc42d8a0b3242182f,0x861d86002cac84a3,0x55d3584fa2c0c09c,0x31e18ae1d869d345},
{0xd402985a6aebab7c,0xe2484e8b9d798d68,0x3ff9c1e7f10c1c5d,0x4d9d50ef1a188d38},
{0xf6f647a6bd61d621,0xd51c337862598d05,0x9ed8fc6d1186caaf,0x5941946bf94019c0},
{0x13fbefc2a3b0e291,0x52d4c35c45380a2d,0x923cfc9cc8e0fed6,0x2c8b7ba8d9f874a8},
{0x9c2f4b81bc9a9ac8,0x372e20bc608968df,0xf7d39da955f3ea72,0x9bc350482e2c686a},
{0xa16c960e23e9037e,0x51922b7a6d39fa7d,0xf0d079380f4a921b,0x968fadae1c073f0c},
{0xab881ac0f1070307,0x26ac16d7f381096a,0xa271eb46db5227b2,0x136478ed3f8133ae},
{0xa54227d10ef8fcf8,0xab03ee2156d8495e,0xfea8b5d5c206789,0x8fa49c7876f775e2},
{0xaf5dac83dc16fc81,0x801dd97edd1f584b,0xc18bfd6c2827fd20,0xc7967b79a716a83},
{0xb49af71043656537,0x9a81e43ce9cfe9e9,0xba88d8fae17ea4c9,0x745c51d884c4125},
{0x3cce52cf5c4f1d6e,0x7edb419d0521489c,0x201f7a076e919065,0x767d99bcdc8034e8},
{0x59d3faeb429e29de,0xfc93d180e7ffc5c3,0x13837a3725ebc48b,0x49c780f9bd388fd0},
{0x7cc7aa3795145483,0xef67b66dacdfc560,0x7262b4bc466672dd,0x556bc4769c601c58},
{0x8c9cb886cdbde7d0,0x4b927ef91dacce25,0x5c891e5494b1ce9f,0x71278a83de0ed64b},
{0x0,0x0,0x0,0x0},
{0xa1948523fffffffe,0xa36009f294b2a592,0x64b8ed486ee51e77,0x46122acb6cf15321},
{0xd402985a6aebab7c,0xe2484e8b9d798d68,0x3ff9c1e7f10c1c5d,0x4d9d50ef1a188d38},
{0x8c9cb886cdbde7d0,0x4b927ef91dacce25,0x5c891e5494b1ce9f,0x71278a83de0ed64b},
{0x13fbefc2a3b0e291,0x52d4c35c45380a2d,0x923cfc9cc8e0fed6,0x2c8b7ba8d9f874a8},
{0x59d3faeb429e29de,0xfc93d180e7ffc5c3,0x13837a3725ebc48b,0x49c780f9bd388fd0},
{0xa16c960e23e9037e,0x51922b7a6d39fa7d,0xf0d079380f4a921b,0x968fadae1c073f0c},
{0xb49af71043656537,0x9a81e43ce9cfe9e9,0xba88d8fae17ea4c9,0x745c51d884c4125},
{0xa54227d10ef8fcf8,0xab03ee2156d8495e,0xfea8b5d5c206789,0x8fa49c7876f775e2},
{0xa54227d10ef8fcf8,0xab03ee2156d8495e,0xfea8b5d5c206789,0x8fa49c7876f775e2},
{0xb49af71043656537,0x9a81e43ce9cfe9e9,0xba88d8fae17ea4c9,0x745c51d884c4125},
{0xa16c960e23e9037e,0x51922b7a6d39fa7d,0xf0d079380f4a921b,0x968fadae1c073f0c},
{0x59d3faeb429e29de,0xfc93d180e7ffc5c3,0x13837a3725ebc48b,0x49c780f9bd388fd0},
{0x13fbefc2a3b0e291,0x52d4c35c45380a2d,0x923cfc9cc8e0fed6,0x2c8b7ba8d9f874a8},
{0x8c9cb886cdbde7d0,0x4b927ef91dacce25,0x5c891e5494b1ce9f,0x71278a83de0ed64b},
{0xd402985a6aebab7c,0xe2484e8b9d798d68,0x3ff9c1e7f10c1c5d,0x4d9d50ef1a188d38},
{0xa1948523fffffffe,0xa36009f294b2a592,0x64b8ed486ee51e77,0x46122acb6cf15321},
{0x0,0x0,0x0,0x0}
};