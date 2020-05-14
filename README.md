# <center>VIFB：A Visible and Infared Image Fusion Benchmark</center>
This is the official webpage for VIFB. To the best of our knowledge, VIFB is the first (and the only one to date) benchmark in the field of visible-infrared image fusion (VIF), aiming to provide a platform to perform fair and comprehensive performance comparision of VIF algorithms. Currently, **21 image pairs, 20 fusion algorithms and 13 evaluation metrics** are integrated in VIFB, which can be utilize to compare performances. All the fusion results are also available that can be used by users directly.

For more details, please refer to the following paper:

**VIFB：A Visible and Infared Image Fusion Benchmark**  
Xingchen Zhang, Ping Ye, Gang Xiao  
From Shanghai Jiao Tong University & Imperial College London  
Contact: xingchen.zhang@imperial.ac.uk  
[[Download paper](https://arxiv.org/abs/2002.03322)]    

### Abstract
Visible and infrared image fusion is an important area in image processing due to its numerous applications. While much progress has been made in recent years with efforts on developing fusion algorithms, there is a lack of code library and benchmark which can gauge the state-of-the-art. In this paper, after briefly reviewing recent advances of visible and infrared image fusion, we present a visible and infrared image fusion benchmark (VIFB) which consists of **21 image pairs, a code library of 20 fusion algorithms and 13 evaluation metrics**. We also carry out large scale experiments within the benchmark to understand the performance of these algorithms. By analyzing qualitative and quantitative results, we identify effective algorithms for robust image fusion and give some observations on the status and future prospects of this field.

### Methods integrated
1. ADF [1]
2. CBF [2]
3. CNN [3]
4. DLF [4]
5. FPDE [5]
6. GFCE [6]
7. GFF [7]
8. GTF [8]
9. HMSD_GF [6]
10. Hybrid_MSD [9]
11. IFEVIP [10]
12. LatLRR [11]
13. MGFF [12]
14. MST_SR [13]
15. MSVD [14]
16. NSCT_SR [13]
17. ResNet [15]
18. RP_SR [13]
19. TIF [16]
20. VSMWLS [17]


### Evaluation metrics integrated
1. Avgerage gradient
2. Cross entropy
3. Edge intensity
4. Entropy
5. Mutual information
6. PSNR
7. Qabf
8. Qcb
9. Qcv
10. RMSE
11. Spatial frequency
12. SSIM
23. SD 


### Examples of fused images
![](https://github.com/xingchenzhang/Visible-infrared-image-fusion-benchmark/blob/master/fusion-fight.jpg)

## Citation
**If you find this work useful, please consider citing**:
    
    @misc{zhang2020vifb,
    title={VIFB: A Visible and Infrared Image Fusion Benchmark},
    author={Xingchen Zhang and Ping Ye and Gang Xiao},
    year={2020},
    eprint={2002.03322},
    archivePrefix={arXiv},
    primaryClass={cs.CV}}


### References
[1] Durga Prasad Bavirisetti and Ravindra Dhuli. Fusion of infrared and visible sensor images based on anisotropic diffusion and karhunen-loeve transform. IEEE Sensors Journal,
16(1):203–209, 2016.  
[2] B. K. Shreyamsha Kumar. Image fusion based on pixel significance using cross bilateral filter. Signal, Image and Video
Processing, 9(5):1193–1204, Jul 2015  
[3] Yu Liu, Xun Chen, Juan Cheng, Hu Peng, and Zengfu Wang. Infrared and visible image fusion with convolutional neural
networks. International Journal of Wavelets, Multiresolution and Information Processing, 16(03):1850018, 2018  
[4] Hui Li, Xiao-Jun Wu, and Josef Kittler. Infrared and visible image fusion using a deep learning framework. 24th
International Conference on Pattern Recognition, 2018.  
[5] Durga Prasad Bavirisetti, Gang Xiao, and Gang Liu. Multisensor image fusion based on fourth order partial differential equations. In 2017 20th International Conference on
Information Fusion (Fusion), pages 1–9. IEEE, 2017.  
[6] Zhiqiang Zhou, Mingjie Dong, Xiaozhu Xie, and Zhifeng Gao. Fusion of infrared and visible images for night-vision
context enhancement. Applied optics, 55(23):6480–6490, 2016.    
[7] Shutao Li, Xudong Kang, and Jianwen Hu. Image fusion with guided filtering. IEEE Transactions on Image
processing, 22(7):2864–2875, 2013.  
[8] Jiayi Ma, Chen Chen, Chang Li, and Jun Huang. Infrared and visible image fusion via gradient transfer and total variation
minimization. Information Fusion, 31:100–109, 2016.   
[9] Zhiqiang Zhou, Bo Wang, Sun Li, and Mingjie Dong. Perceptual fusion of infrared and visible images through a hybrid multi-scale decomposition with gaussian and bilateral
filters. Information Fusion, 30:15–26, 2016.  
[10] Yu Zhang, Lijia Zhang, Xiangzhi Bai, and Li Zhang. Infrared and visual image fusion through infrared feature extraction
and visual information preservation. Infrared Physics & Technology, 83:227 – 237, 2017.  
[11] Hui Li and Xiaojun Wu. Infrared and visible image fusion using latent low-rank representation. arXiv preprint
arXiv:1804.08992, 2018.  
[12] Durga Prasad Bavirisetti, Gang Xiao, Junhao Zhao, Ravindra Dhuli, and Gang Liu. Multi-scale guided image and video
fusion: A fast and efficient approach. Circuits, Systems, and Signal Processing, 38(12):5576–5605, Dec 2019.   
[13] Yu Liu, Shuping Liu, and Zengfu Wang. A general framework for image fusion based on multi-scale transform and
sparse representation. Information Fusion, 24:147–164, 2015.  
[14] VPS Naidu. Image fusion technique using multi-resolution singular value decomposition. Defence Science Journal,
61(5):479–484, 2011.  
[15] Hui Li, Xiao-Jun Wu, and Tariq S Durrani. Infrared and visible image fusion with resnet and zero-phase component
analysis. Infrared Physics & Technology, 102:103039, 2019.  
[16] Durga Prasad Bavirisetti and Ravindra Dhuli. Two-scale image fusion of visible and infrared images using saliency detection. Infrared Physics & Technology, 76:52–64, 2016.  
[17] Jinlei Ma, Zhiqiang Zhou, Bo Wang, and Hua Zong. Infrared and visible image fusion based on visual saliency map
and weighted least square optimization. Infrared Physics & Technology, 82:8–17, 2017.  