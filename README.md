# SP_Face_Detection
執行方式：
1. 執行skin_filter_data_collect.m，收集training data
2. 執行skin_filter_normalize.m，對收集到的training data做normalization
3. 執行skin_filter_train.m，利用svm訓練模型
4. 執行skin_filter_test.m，測試訓練出的模型並同時增加training data以增加新一輪訓練的精準度，並對得到skin filter結果做後處理
5. 執行skin_filter_binary_seg.m，對偵測結果做二值化分割，呈現結果如skin_filter_final_result所示
6. 執行ellipse_matching.m，排除掉skin filter偵測為膚色，但形狀不是橢圓形的區塊，執行結果如ellipse_matching_result中圖片所示
7. 執行eyemap.m，篩選出圖片中可能是眼睛的區域
8. 執行mouthmap.m，篩選出圖片中可能是嘴巴的區域
9. 執行feature_map_summary.m，綜合位置、顏色等資訊，判斷出一張圖片中是否有人臉存在