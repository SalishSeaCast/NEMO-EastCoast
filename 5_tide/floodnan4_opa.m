function data_new=floodnan4_opa(data, mask, bx)
% data_new=floodnan(data, mask, bx)
%
% data - data field with NaNs over land (and MISSING values)
% mask - true land mask with NaNs over land, zeros over water
% bx   - box size
%
% January 28, 2003
% lukeman@phys.ocean.dal.ca

[si,sj]=size(data);
frame=ones(si+bx*2,sj+bx*2)*NaN;                                % create frame around data
data_framed=frame; data_framed(bx+1:bx+si,bx+1:bx+sj)=data;     % and mask so that cellbox
mask_framed=frame; mask_framed(bx+1:bx+si,bx+1:bx+sj)=mask;     % indices always in range

%[ii jj]=size(data_framed)

% lat=ones(si+bx*2,sj+bx*2)*NaN;
% lat(bx+1:bx+si,bx+1:bx+sj)=lata;
%lat=[ones(1,bx)*NaN lat' ones(1,bx)*NaN];

[x,y]=meshgrid(-bx:bx,-bx:bx);

Idata=find(isnan(data_framed));                                 % get index of data NaNs
if isempty(Idata),
    error('land mask on data must be NaNs')
end

Imask=find(isnan(mask_framed));                                 % get index of mask NaNs
if isempty(Imask),
    error('land mask on mask must be NaNs')
end

Idiff1=setdiff(Imask,Idata);                % find new land points
data_framed(Idiff1)=NaN;                    % and set as land

Idiff2=setdiff(Idata,Imask);                % find points to be flooded

data_new=data_framed;

Idiff_again=[];

for again=1:2;
    
%again=1;
%while again,
    count=0;
for n=1:length(Idiff2),
    n;
    %disp([ 'n= ' num2str(n) 'Idiff2' num2str(length(Idiff2))])
    [i,j] = ind2sub(size(data_framed),Idiff2(n));               % get i,j index from linear one
    %disp(['i' num2str(i)  'j' num2str(j)])
    cellbox=data_framed(i-bx:i+bx,j-bx:j+bx);                   % extract box around point
    cbnan=~isnan(cellbox);
    
    
    if min(isnan(cellbox(:))),                                  % if box is all NaNs then
        count=count+1;                                          % keep track of index as
        Idiff_again(count)=Idiff2(n);                           % we'll have to pass again
    else
        %R=abs(x.*cos(lat(i,j)*pi/180) + y.*sqrt(-1));
        R=abs(x + y.*sqrt(-1));
        weights=exp(-R./.5);
        
  % size(weights)
   %size(R)
   
        data_new(i,j)=sum(cellbox(cbnan).*weights(cbnan))./sum(weights(cbnan)); % take weighted average    
       
    end
end
    if length(Idiff_again>=1),                  % set up to do for loop again for cellboxes that              
        Idiff2=Idiff_again;                     % were all NaN in previous iteration.
        Idiff_again=[];
        data_framed=data_new;
    else
        %again=0;                                % otherwise, we're done.
    end
end


data_new=data_new(bx+1:end-bx,bx+1:end-bx);     % remove NaN frame from output
