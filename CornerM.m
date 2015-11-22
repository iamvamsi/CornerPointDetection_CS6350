% ========================================================+
% Authors:                                                |
% --------------------------------------------------------+
% Vamsi Krishan D (dvamsi@cse.iitm.ac.in)    
% Venkatehs G (gvenky@cse.iitm.ac.in)
% ========================================================+
% File : CornerDeector.m                                  |
% ========================================================+
% Corner Detection :                                      |
% --------------------------------------------------------+
% AN EFFECTIVE APPROACH TO CORNER POINT DETECTION THROUGH
% MULTIRESOLUTION ANALYSIS -- By Yang Bai, Hairong Qi     |
% --------------------------------------------------------+

clear;
input_RGB=imread('img4.ppm');
%Particular image to load
input_gray=rgb2gray(input_RGB);
%Converting to gray scale
S=size(input_gray);
%PG_size is the size of the Polarized Gaussian kernel used.
if S(1)<S(2)
    PG_size=10
else
    PG_size=10
end
%%%%Generating Binary masks
%%%%Using analytical properties to generate Binary masks.
P1 = true(floor(PG_size/2),PG_size);
P1_add=(-1)*true(PG_size-floor((PG_size/2)),PG_size);
P1=cat(1,P1,P1_add);
id1=floor(PG_size/sqrt(3));
id2=PG_size;
P2=-1*ones(id1,id2);
for i=2:PG_size
    j=1;
   while((id1-j)/i) > (1/sqrt(3))
       P2(j,i)=1;
       j=j+1;
   end
end
P2(find(P2(:,1)==-1))=1;
P2=cat(1,P2,-1*ones(floor((PG_size-id1)/2),PG_size));
P2=cat(1,ones(PG_size-floor((PG_size-id1)/2)-id1,PG_size),P2);
id1=floor(PG_size/sqrt(3))
id2=PG_size;
P3=-1*ones(id2,id1);
for i=2:id1
    j=1;
   while((id2-j)/i)> (sqrt(3))
       P3(j,i)=1;
       j=j+1;
   end
end
P3=cat(2,P3,-1*ones(PG_size,floor((PG_size-id1)/2)));
P3=cat(2,ones(PG_size,PG_size-floor((PG_size-id1)/2)-id1),P3);
P4=true(PG_size,floor(PG_size/2));
P4_add=(-1)*true(PG_size,PG_size-floor((PG_size/2)));
P4=horzcat(P4,P4_add);
P5=flipud(P3);
P6=flipud(P2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%Generating polarized gaussian kernels
PG1=P1.*fspecial('gaussian',PG_size);
PG2=P2.*fspecial('gaussian',PG_size);
PG3=P3.*fspecial('gaussian',PG_size);
PG4=P4.*fspecial('gaussian',PG_size);
PG5=P5.*fspecial('gaussian',PG_size);
PG6=P6.*fspecial('gaussian',PG_size);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Res(S(1),S(2))=0;
for N=1:4
    [Ci,Si] = wavedec2(input_gray,N,'haar');
    
    [chd,cvd,cdd] = detcoef2('all',Ci,Si,N);
    temp=(abs(chd.*cvd.*cdd)).^(1/3);
    for it=1:N
        temp=impyramid(temp,'expand');
    end;
    Itemp=imresize(temp,S,'bilinear');
    Res = Res + Itemp;
    clear temp;
end   
 C_m=imregionalmax(Res);
 %%%C_m contains both true and false corners.
 C_m=double(C_m);
 [row,col]= find(C_m==1);
 C_con = conv2(C_m,fspecial('gaussian',7));
for i= 1:size(row)
    C_m(row(i),col(i))=1./C_con(row(i),col(i));  %%%Approximation to corner scale map , which gives less weightage to false corners (which are generally crowded)
end
%%%%%%%%%%%Using polarized gaussian kernels to achieve rotational invariance.
%%%%%%%%%%%To separate edge points and some real corners
%%%% with a smaller scale parameter, convolve these points
%%%% with the 6 polarized Gaussian kernels, respectively
R1=conv2(C_m,PG1,'same');
R2=conv2(C_m,PG2,'same');
R3=conv2(C_m,PG3,'same');
R4=conv2(C_m,PG4,'same');
R5=conv2(C_m,PG5,'same');
R6=conv2(C_m,PG6,'same');
%%%%Remove corners that have small responses for all 6 polarized
%%%%Gaussian kernels as these would indicate the edge points.
THRESHOLD=0.97*(min(min(R1))+min(min(R2))+min(min(R3))+min(min(R4))+min(min(R5))+min(min(R6)));
j=1;
for i=1:size(row)
    Th = R1(row(i),col(i))+R2(row(i),col(i))+R3(row(i),col(i))+R4(row(i),col(i))+R5(row(i),col(i))+R6(row(i),col(i));
    if(Th >= THRESHOLD)
        Nrow(j)=row(i);
        Ncol(j)=col(i);
        j=j+1;
    end
end
%%%%%%Displaying them on the original image
imshow(input_RGB);
hold on;
scatter(Nrow,Ncol,3,'filled');
hold off;