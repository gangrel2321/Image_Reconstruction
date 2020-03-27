%function ReadData

I1 = double(imread('barbara_contaminated.png'));
I2 = double(imread('cameraman_contaminated.png'));
ground_camera = double(imread('cameraman.png'));
image = I1; %specify which image you want to use

ground_cam = ground_camera(:);

figure;
imshow(image,[]);
load Omega  %%%%% I1(Omega), I2(Omega) are not contaminated

%%%%%%   Wavelet transformation  = W*u, W^T*W = I
%%%%%  Inverse wavelet transformation  = W^T*u

coef = swt2(image,1,1);
total = 1:256*256;
Omega_C = setdiff(total, Omega);
omega_size = size(Omega,1);

first_layer = coef(:,:,1); %used for inverse wavelet, don't threshold
lambda_k = zeros(256*256*8,1);
theta_k = zeros(size(Omega,1),1);
x_k = image(:);
q_k = coef_to_vec(coef,2);

trans_x_k = reshape(x_k, [], 256); %make x_k square pixel matrix again
wt = swt2(trans_x_k,1,1); %calculate wavelet transform
wt = coef_to_vec(wt,2); 

mu_1 = 150;
mu_2 = 150;
epsilon = 17;

%create a sparse A matrix
A = sparse(omega_size,256*256); 
for i=1:omega_size
    A(i,Omega(i)) = 1; 
end
b = x_k(Omega);

update_term = zeros(256,256,9);

for i=1:300
   x_old = x_k;
   q_old = q_k;
   theta_old = theta_k;

   %argmin_x, solve for x_k
   update_term(:,:,2:end) = (reshape(q_k, [256,256,8]) - reshape(lambda_k, [256,256,8]));
   update_term(:,:,1) = first_layer; 
   w_trans_term = iswt2( update_term ,1,1); %w' * (mu*q_k - lambda)
   w_trans_term = w_trans_term(:); 
   x_k = conjugate((mu_1*(A'*A) + mu_2*speye(size(A',1)) ), (mu_1*A'*(b - theta_k) + mu_2*w_trans_term));

   %calculate wavelet transform of new x_k
   trans_x_k = reshape(x_k, [], 256); %make x_k square pixel matrix again
   wt = swt2(trans_x_k,1,1); %calculate wavelet transform
   wt = coef_to_vec(wt,2); 
   
   %soft thresholding
   q_k = sign(wt - (lambda_k/mu_2)) .* max(abs(wt - (lambda_k/mu_2)) - 1, 0); 
   
   % use l2 norm as stopping condition
   diff = norm(q_k - q_old,2);
   if diff < epsilon
       fprintf('Result Met\n')
       i
       break
   end
   %update theta and lambda
   
   theta_k = theta_k + (x_k(Omega) - b); %Ax_k = x_k(Omega)
   lambda_k = lambda_k + (wt - q_k);
    
   %print out progress (for testing)
   cur_Im = trans_x_k;
   if mod(i,10) == 0
       figure;
       imshow(cur_Im, []);
   end
   
   %q = norm(q_k,2)
   %x = norm(x_k - x_old,2)
   %lambda = norm(lambda_k, 2)
   %theta = norm(theta_k, 2)
   %wt_norm = norm(wt, 2)

end
close all; 
figure; 
imshowpair(image, cur_Im, 'montage')
