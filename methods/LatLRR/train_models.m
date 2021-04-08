% train different project matrices with different size 8,16 
% training data in 'training_data'
training_path = './training_data/';
col = 1200;
unit = 8;
L_8 = latent_matrix(training_path, col, unit);
save('./models/L_8.mat', 'L_8');
unit = 16;
L_16 = latent_matrix(training_path, col, unit);
save('./models/L_16.mat', 'L_16');




