# <center>VIFB：A Visible and Infared Image Fusion Benchmark</center>
This is the official webpage of VIFB. 

**VIFB is the first (and the only one to date) benchmark in the field of visible-infrared image fusion (VIF)**, aiming to provide a platform to perform fair and comprehensive performance comparision of VIF algorithms. Currently, **21 image pairs, 20 fusion algorithms and 13 evaluation metrics** are integrated in VIFB, which can be utilized to compare performances conveniently. All the fusion results are also available that can be used by users directly. In addition, more test images, fusion algorithms (in Matlab), evaluation metrics and fused images can be easily added using the provided toolkit.

For more details, please refer to the following paper:

**VIFB: A Visible and Infared Image Fusion Benchmark**  
Xingchen Zhang, Ping Ye, Gang Xiao  
In the Proceedings of CVPR Workshop, 2020  
From Shanghai Jiao Tong University & Imperial College London  
Contact: xingchen.zhang@imperial.ac.uk  
[[Download paper](https://arxiv.org/abs/2002.03322)]

Chinese readers can also refer to [[this link](https://mp.weixin.qq.com/s/KB-f8maHuWZLUbvbxUrPxw)] for more details and the motivation of this benchmark. We really hope this benchmark can contribute to the development of image fusion field. 

**If you use this code, please cite**:

	@inproceedings{zhang2020vifb,
	title={VIFB: A Visible and Infrared Image Fusion Benchmark},
	author={Zhang, Xingchen and Ye, Ping and Xiao, Gang},
	booktitle={Proceedings of the IEEE/CVF Conference on Computer Vision and Pattern Recognition Workshops},
	year={2020}}  

	@article{zhang2023visible,
	title={Visible and Infrared Image Fusion Using Deep Learning},
	author={Zhang, Xingchen and Demiris, Yiannis},
	journal={IEEE Transactions on Pattern Analysis and Machine Intelligence},
	year={2023},
	publisher={IEEE}}

## Abstract
Visible and infrared image fusion is an important area in image processing due to its numerous applications. While much progress has been made in recent years with efforts on developing fusion algorithms, there is a lack of code library and benchmark which can gauge the state-of-the-art. In this paper, after briefly reviewing recent advances of visible and infrared image fusion, we present a visible and infrared image fusion benchmark (VIFB) which consists of **21 image pairs, a code library of 20 fusion algorithms and 13 evaluation metrics**. We also carry out large scale experiments within the benchmark to understand the performance of these algorithms. By analyzing qualitative and quantitative results, we identify effective algorithms for robust image fusion and give some observations on the status and future prospects of this field.

## What are contained in VIFB

### Dataset
The dataset in VIFB is a test set, which is collected by the authors from the [Internet](https://www.ino.ca/en/solutions/video-analytics-dataset/) and from fusion tracking datasets [1,2,3]. Image registration between RGB and infrared images are also conducted by corresponding authors. We appreciate the authors of these datasets very much for making these images publicly available for research. **Please also cite these papers and the [link](https://www.ino.ca/en/solutions/video-analytics-dataset/) if you use VIFB**. Thanks!

![](https://github.com/xingchenzhang/Visible-infrared-image-fusion-benchmark/blob/master/dataset.jpg)


### Methods integrated
Currently, we have integrated 20 VIF algorithms in VIFB. Many thanks for the authors of these algorithms for making their codes available to the community. **Please cite these papers as well if you use VIFB**. Thanks!

1. ADF [4]
2. CBF [5]
3. CNN [6]
4. DLF [7]
5. FPDE [8]
6. GFCE [9]
7. GFF [10]
8. GTF [11]
9. HMSD_GF [9]
10. Hybrid_MSD [12]
11. IFEVIP [13]
12. LatLRR [14]
13. LP_SR [15]
14. MGFF [16]
15. MSVD [17]
16. NSCT_SR [15]
17. ResNet [18]
18. RP_SR [15]
19. TIF [19]
20. VSMWLS [20]

The download links of each algorithm can be found on [this website](https://zhuanlan.zhihu.com/p/342971809). For each algorithm, we use original settings reported by corresponding authors in their papers. For deep learning-based methods, the pretrained model provided by corresponding authors are used. We did not train these algorithms.

We modified the interfaces of these algorithms to run them in VIFB conveniently. For other methods written in MATLAB, they can also be added to VIFB by chaning the interface. For algorithms written in Python or other languages, we suggest the users change the name of the fused images and put them in the output folder. Then the evaluation methods can be computed by finishing correponding settings.

#### Updated: we added results for 5 VIF methods:
1. IFCNN [21]
2. SeAFusion [22]
3. SwinFusion [23]
4. U2Fusion [24]
5. YDTR [25]

### Evaluation metrics integrated
We have integrated 13 evaluation metrics in VIFB. The codes are collected from the Internet, forum, etc. and checked by the authors.

Many thanks to the authors of these evaluation metric codes for sharing their codes with the community. Please forgive us because we can not find the original source of these codes so we can not mention everyone here. Please let us know if you think you are the authors of these evaluation metric codes and we will add references and acknowledgements. Many thanks!

1. Average gradient
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
13. SD 


### Examples of fused images
![](https://github.com/xingchenzhang/Visible-infrared-image-fusion-benchmark/blob/master/fusion-fight.jpg)


## How to use
### How to run algorithms
1. Please set the algorithms you want to run in util\configMethods.m
2. Please set the images you want to fuse in util\configImgsVI.m and util/configImgsIR.m, and change the path of these images
3. DLF and ResNet methods need MatConvNet to run. One should set the path to MatConvNet in run_ DLF.m (line 28) nad run_ ResNet.m (line 10), respectively.
4. DLF requires 'imagenet-vgg-verydeep-19.mat' to run. Please download it and put it inside methods\DLF
5. ResNet requires 'imagenet-resnet-50-dag.mat' to run. Please download it and put it inside methods\ResNet\models
6. To run GFF, please set your own path in run_GFF.m (line 17) 
7. To run LP_SR for the first time, please run "make" in the path "...\methods\LP_SR\sparsefusion\ksvdbox\ompbox\private". Similarly, to run RP_SR for the first time, run "make" in the path "...\methods\RP_SR\sparsefusion\ksvdbox\ompbox\private". To run NSCT_SR for the first time, run "make" in the path "...\methods\NSCT_SR\sparsefusion\ksvdbox\ompbox\private".
8. main_running.m is used to run the fusion algorithms. Please change the output path in main_running.m.
9. Enjoy!


### How to compute evaluation metrics
1. Please set the metrics you want to compute in util\configMetrics.m
2. compute_metrics.m is used to compute evaluation metrics. Please change the output path in compute_metrics.m
3. Enjoy!

### How to add algorithms (or fused images)
1. For methods written in MATLAB, please put them in the folder methods. For example, for method "ADF", put the codes inside a folder called "ADF", and put the folder "ADF" inside "Methods". Then change the main file of ADF to run_ADF.m. In run_ADF.m, please change the interface as according to examples provided in VIFB.
2. For algorithms written in Python or other languages, we suggest the users change the name of the fused images  according to examples provided and put them in the output folder. Then add the methods in util\configMethods.m. Then, the evaluation metrics can be computed.

## Acknowledgement
The overall framework of VIFB is created based on OTB [26]. We thank the authors of OTB very much for making OTB publicly available. We also thank all authors of the integrated images, VIF methods and evaluation metrics for sharing their work to the community! 

### References
[1] C. Li, X. Liang, Y. Lu, N. Zhao, and J. Tang, “Rgb-t object tracking: benchmark and baseline,” Pattern Recognition, p106977, 2019.  
[2] J. W. Davis and V. Sharma, “Background-subtraction using contour-based fusion of thermal and visible imagery,”
Computer vision and image understanding, vol. 106, no. 2-3, pp. 162–182, 2007.  
[3] C. O’Conaire, N. E. O’Connor, E. Cooke, and A. F. Smeaton, “Comparison of fusion methods for thermo-visual surveillance tracking,” in 2006 9th International Conference on Information Fusion. IEEE, 2006, pp. 1–7.  
[4] Durga Prasad Bavirisetti and Ravindra Dhuli. Fusion of infrared and visible sensor images based on anisotropic diffusion and karhunen-loeve transform. IEEE Sensors Journal,
16(1):203–209, 2016.      
[5] B. K. Shreyamsha Kumar. Image fusion based on pixel significance using cross bilateral filter. Signal, Image and Video
Processing, 9(5):1193–1204, Jul 2015  
[6] Yu Liu, Xun Chen, Juan Cheng, Hu Peng, and Zengfu Wang. Infrared and visible image fusion with convolutional neural
networks. International Journal of Wavelets, Multiresolution and Information Processing, 16(03):1850018, 2018  
[7] Hui Li, Xiao-Jun Wu, and Josef Kittler. Infrared and visible image fusion using a deep learning framework. 24th
International Conference on Pattern Recognition, 2018.  
[8] Durga Prasad Bavirisetti, Gang Xiao, and Gang Liu. Multisensor image fusion based on fourth order partial differential equations. In 2017 20th International Conference on
Information Fusion (Fusion), pages 1–9. IEEE, 2017.  
[9] Zhiqiang Zhou, Mingjie Dong, Xiaozhu Xie, and Zhifeng Gao. Fusion of infrared and visible images for night-vision
context enhancement. Applied optics, 55(23):6480–6490, 2016.    
[10] Shutao Li, Xudong Kang, and Jianwen Hu. Image fusion with guided filtering. IEEE Transactions on Image
processing, 22(7):2864–2875, 2013.  
[11] Jiayi Ma, Chen Chen, Chang Li, and Jun Huang. Infrared and visible image fusion via gradient transfer and total variation
minimization. Information Fusion, 31:100–109, 2016.   
[12] Zhiqiang Zhou, Bo Wang, Sun Li, and Mingjie Dong. Perceptual fusion of infrared and visible images through a hybrid multi-scale decomposition with gaussian and bilateral
filters. Information Fusion, 30:15–26, 2016.  
[13] Yu Zhang, Lijia Zhang, Xiangzhi Bai, and Li Zhang. Infrared and visual image fusion through infrared feature extraction
and visual information preservation. Infrared Physics & Technology, 83:227 – 237, 2017.  
[14] Hui Li and Xiaojun Wu. Infrared and visible image fusion using latent low-rank representation. arXiv preprint
arXiv:1804.08992, 2018.  
[15] Yu Liu, Shuping Liu, and Zengfu Wang. A general framework for image fusion based on multi-scale transform and
sparse representation. Information Fusion, 24:147–164, 2015.   
[16] Durga Prasad Bavirisetti, Gang Xiao, Junhao Zhao, Ravindra Dhuli, and Gang Liu. Multi-scale guided image and video fusion: A fast and efficient approach. Circuits, Systems, and Signal Processing, 38(12):5576–5605, Dec 2019.   
[17] VPS Naidu. Image fusion technique using multi-resolution singular value decomposition. Defence Science Journal,
61(5):479–484, 2011.  
[18] Hui Li, Xiao-Jun Wu, and Tariq S Durrani. Infrared and visible image fusion with resnet and zero-phase component
analysis. Infrared Physics & Technology, 102:103039, 2019.  
[19] Durga Prasad Bavirisetti and Ravindra Dhuli. Two-scale image fusion of visible and infrared images using saliency detection. Infrared Physics & Technology, 76:52–64, 2016.  
[20] Jinlei Ma, Zhiqiang Zhou, Bo Wang, and Hua Zong. Infrared and visible image fusion based on visual saliency map and weighted least square optimization. Infrared Physics & Technology, 82:8–17, 2017.  
[21] Zhang, Y., Liu, Y., Sun, P., Yan, H., Zhao, X., & Zhang, L.. IFCNN: A general image fusion framework based on convolutional neural network. Information Fusion, 54, 99-118, 2020.  
[22] L. Tang, J. Yuan, & J. Ma. Image fusion in the loop of high-level vision tasks: A semantic-aware real-time infrared and visible image fusion network. Information Fusion, 82, 28-42, 2022.  
[23] J. Ma, L. Tang, F. Fan, J. Huang, X. Mei, & Y. Ma. SwinFusion: Cross-domain long-range learning for general image fusion via swin transformer. IEEE/CAA Journal of Automatica Sinica, 9(7), 1200-1217, 2022.  
[24] H. Xu, J. Ma, J. Jiang, X. Guo, & H. Ling. U2Fusion: A unified unsupervised image fusion network. IEEE Transactions on Pattern Analysis and Machine Intelligence, 44(1), 502-518, 2022.  
[25] W. Tang, F. He, & Y. Liu. YDTR: infrared and visible image fusion via y-shape dynamic transformer. IEEE Transactions on Multimedia, 2022.  
[26] Y. Wu, J. Lim, & M. H., Yang. Online object tracking: A benchmark. In Proceedings of the IEEE conference on computer vision and pattern recognition, pp. 2411-2418, 2013.