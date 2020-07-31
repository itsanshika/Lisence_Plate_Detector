% Reading Input Image
A = imread('C:\Users\anshi\Desktop\car2.jpg');
A=imnoise(A,'salt & pepper',0.01);
% RGB to GRAYSCALE CONVERSION
R=A(:, :, 1);
G=A(:, :, 2);
B=A(:, :, 3);
[M, N, ~]=size(A);
gray_img=zeros(M, N, 'uint8');  %matrix to image conversion mat2gray
for i=1:M
         for j=1:N
               gray_img(i, j)=(R(i, j)*0.2989)+(G(i, j)*0.5870)+(B(i, j)*0.114);  %best  possible combination
              end
end
     
  figure(1), imshow(gray_img);
  D=imresize(gray_img,[500,800]);
 
  %Median filtering
output = zeros(500,800);
output = uint8(output);
for i = 1:500
    for j = 1:800  %intesity of pixel in the noisy image is given as noisy(i,j)
        % here we define max and minimum values x and y coordinates of any
        % pixel can take
        xmin = max(1,i-1); % minimum x coordinate has to be greater than or equal to 1
        xmax = min(500,i+1);
        ymin = max(1,j-1);
        ymax = min(800,j+1);
        % the neighbourhood matrix will then be
        temp = D(xmin:xmax, ymin:ymax);
        %now the new intensity of pixel at (i,j) will be median of this
        %matrix
        output(i,j) = median(temp(:));
    end
end
figure(7), imshow(output);
 

% BINARIZATION OF IMAGE
U=output;
[x,y,z]=size(output);

 U=double(U);
     
  sum=0;
    for i=1:x
        for j=1:y
            sum=sum+U(i, j);
        end
    end
   
    level=sum/(x*y);
    binary=zeros(x, y);
   
  for i=1:x
      for j=1:y
        if U(i, j)>= level
                binary(i, j) =1;
        else
            binary(i, j)=0;
        end
     end
  end
   
    figure(2), imshow(binary);
   
 



   %DILATION & EROSION OF IMAGE
   I=binary;
   
   % create structuring element              
se=ones(3,3);  
 
% store number of rows in P and number of columns in Q.            
[P, Q]=size(se);  
 
% create a zero matrix of size I.        
In=zeros(size(I, 1), size(I, 2));
In1=zeros(size(I, 1), size(I, 2));
 
for i=ceil(P/2):size(I, 1)-floor(P/2)
    for j=ceil(Q/2):size(I, 2)-floor(Q/2)
 
        % take all the neighbourhoods.
        on=I(i-floor(P/2):i+floor(P/2), j-floor(Q/2):j+floor(Q/2));  
         
        % take logical se
        nh=on(logical(se));    
 
        % compare and take minimum value of the neighbor  
        % and set the pixel value to that minimum value.    
        In(i, j)=max(nh(:,:));
        In1(i, j)=min(nh(:,:));
    end
end
  figure(3), imshow(In);
  figure(4), imshow(In1);
  In2=In-In1;
  figure(5), imshow(In2);
 
 
    %edge enhancement
kernel = -1*ones(3);
kernel(2,2) = 10;
I = eeconv2(kernel,In2);
figure(8), imshow(I);
 
 
   function result = eeconv2(kernel,targetMat)
%----------------------------------
%Part One: create the padded matrix
%----------------------------------
%determine the dimensions of the two matrices
[m1,n1] = size(kernel);
[m2,n2] = size(targetMat);
%determine how much padding is necessary
extram = (m1 - 1); %this is how many rows will need to be added.
extran = (n1 - 1); %this is how many columns will need to be added.
%initialize the padded matrix
padTarget = zeros(m2+extram,n2+extran);
%fill the middle of padTarget with targetMat
padTarget(1+extram/2:extram/2+m2,1+extran/2:extran/2+n2) = targetMat;
%fill the repetitions on each side in
%------------------------------------
for i = 1:extram/2
       
        %top rows
        padTarget(i,1+extran/2:extran/2+n2) = padTarget(1+extram/2,1+extran/2:extran/2+n2);
        %bottom rows
        padTarget(i+extram/2+m2,1+extran/2:extran/2+n2) = padTarget(extram/2+m2,1+extran/2:extran/2+n2);
        %left columns
        padTarget(1+extram/2:extram/2+m2,i) = padTarget(1+extram/2:extram/2+m2,1+extran/2);
        %right columns
        padTarget(1+extram/2:extram/2+m2,i+extran/2+n2) = padTarget(1+extram/2:extram/2+m2,extran/2+n2);
end
%fill the corners of the padded matrix
%-------------------------------------
%topleft
padTarget(1:extram/2,1:extran/2)=padTarget(extram/2,1+extran/2);
%bottomleft
padTarget(extram/2+m2+1:extram+m2,1:extran/2)=padTarget(m2+extram/2,1+extran/2);
%topright
padTarget(1:extram/2,1+n2+extran/2:n2+extran)=padTarget(1+extram/2,n2+extran/2);
%bottomright
padTarget(extram/2+m2+1:extram+m2,1+n2+extran/2:n2+extran)=padTarget(m2+extram/2,n2+extran/2);
%---------------------------------
%Part Two: perform the convolution
%---------------------------------
%flip the rows and columns of the kernel
%---------------------------------------
flipconvMat = zeros(m1,n1); %preallocate flipconvMat
%flip
for i = 1:m1
    for j = 1:n1
        flipconvMat(i,j) = kernel(1+m1-i,1+n1-j);
    end
end
%now element by element perform the convolution
%outer loop selects which pixel in padTarget we are centering on at the moment
result = zeros(m2,n2); %preallocate result
for i = 1+extram/2:extram/2+m2
    for j = 1+extran/2:extran/2+n2
   
        convSum = 0; %the running total
        for m = 1:m1
            for n = 1:n1
                convSum = convSum+flipconvMat(m,n)*padTarget(i-extram/2-1+m,j-extran/2-1+n);
            end
        end
        result(i-extram/2,j-extran/2)=convSum;
    end
end
end