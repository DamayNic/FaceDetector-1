function sc = ApplyDetector(Cparams,ii_im,mu,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%This function compute the sccore of the detector fot the current image.
%%Cparams:      Strong classifier parameters.
%%.             Cparams.alphas: (T,1) 
%%.             Cparams.Thetas: (T,3) -> (index,theta,p)
%%.             Cparams.fmat:  (HW,d)
%%.             Cparams.all_filters: (d,5);
%%ii_im: Integral Image (N,HW)
%%mu: mean of ....  (1,N)   
%%sigma: standart deviation of... (1,N)    

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3 %The integral images have been computed on normalized images.
    True_applied_filters = (Cparams.fmat(:,Cparams.Thetas(:,1))'*ii_im'); 
else  %%The integral images have been not computed on normalized images. The filter have to be changed. 
      %F1*=F1/sigma
      %F2*=F2/sigma
      %F3*=F3/sigma+mu/sigma*h*w
      %F4*=F4/sigma
    
    T=size(Cparams.Thetas,1);
    N=size(ii_im,1);
    
    chosen_filters=Cparams.Thetas(:,1); %(T,1)
    
    Applied_Chosen_filters=Cparams.fmat(:,chosen_filters)'*ii_im'; %(T,N);
    
    IsType3=(Cparams.all_filters(Cparams.Thetas(:,1),5)==3);
    w=Cparams.all_filters(Cparams.Thetas(:,1),3);
    h=Cparams.all_filters(Cparams.Thetas(:,1),4);
    
    True_applied_filters = (Applied_Chosen_filters./repmat(sigma',T,1))+repmat(IsType3.*w.*h,1,N).*repmat((mu./sigma)',T,1);
end

sc = Cparams.alphas'*((repmat(Cparams.Thetas(:,3),1,size(ii_im,1)).*True_applied_filters)<repmat(Cparams.Thetas(:,3).*Cparams.Thetas(:,2),1,size(ii_im,1)));

end

