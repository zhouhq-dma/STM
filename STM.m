function STM(z) % z is the coordination along z of the STM image plane. z is in [0 az]. unit: Ang
ispin=1;    %cartesian positions
repeat=[1 1]; % extend the presentation of charge plot

%%%%%%%%% read CHGCAR %%%%%%%%%%%%%%%
    fid=fopen('PARCHG');
    textscan(fid,'%s',1);
    tmp=textscan(fid,'%f',1);scale=tmp;
	tmp =  textscan(fid,'%f',9);
    lattice=reshape(tmp{1},3,3)';
    tmp=textscan(fid,'',1); nspec=length(tmp);
    disp(['There are ',num2str(nspec),' species.']);
    az=lattice(3,3);
    disp(['Lattice constant along z is :' ,num2str(az)]);
    
     if (isnumeric(tmp{1}))&&(~isempty(tmp{1}))
        num(1:nspec)=tmp{1:nspec};    
    elseif(ischar(tmp{1}))
        type(1:nspec)=tmp{1:nspec};
        disp(['         ', type(1:nspec)]);
        tmp=textscan(fid,'%d',nspec);
        num(1:nspec)=tmp{1:nspec};
    elseif(isempty(tmp{1}))
        tmp=textscan(fid,'%s',nspec); % it is only one field     
        type=tmp{1}; 
        disp(type');
        tmp=textscan(fid,'%d',nspec);
        num(1:nspec)=tmp{1};  
     end    
    nat=sum(num);
    disp(['        ', num2str(num)]);
    disp(['Total number of atoms: ',num2str(nat)]);
   
    tmp=textscan(fid,'%s',1); cartORdirect=tmp{1};
    tmp=textscan(fid,'%f',nat*3); pos=reshape(tmp{1},3,nat)';
 
    tmp= textscan(fid,'%f',3); 
    NGX=tmp{1}(1); NGY=tmp{1}(2); NGZ=tmp{1}(3); ngrid=NGX*NGY*NGZ;
    
    tmp=textscan(fid,'%f',ngrid); 
    chg=reshape(tmp{1}, NGX,NGY,NGZ);
	fclose(fid);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if( nargin ==0) % without input of nz
            aveZ=  sum(pos(:,3))/nat;
            tmp=lattice(3,3);    
            dz=lattice(3,3)/NGZ;
            nz=  floor(21/dz);% floor(( aveZ - min(pos(:,3)))/dz);
            if (nz==0) nz=1; end
    else          
        dz=az/NGZ;
        nz=  floor(mod(z,az)/dz)+1;
        disp(['Input z is ' ,num2str(z)]);
   
    end

    zreal=(nz-1)*dz;
    disp(['Will show slice at z = ' ,num2str(zreal),' Angstrom']);
    
    chgSlice=chg(:,:,nz);
   % draw %
   figure;
%    n1=linspace(0,1,NGX); a1=lattice(1,:);
%    n2=linspace(0,1,NGY); a2=lattice(2,:);
%    for nx=1:NGX
%        for ny=1:NGY
%            xx(nx,ny)= a1(1) * n1(nx) + a2(1) * n2(ny);
%            yy(nx,ny)= a1(2) * n1(nx) + a2(2) * n2(ny);
%        end
%    end
   % one cell
%    surf(xx,yy,-chgSlice); shading interp; view(0,90);
%    axis tight; box on; colormap hot; colorbar

 % extended cell
   n1=linspace(0,1,NGX*repeat(1)); a1=lattice(1,:);
   n2=linspace(0,1,NGY*repeat(2)); a2=lattice(2,:);
   for nx=1:NGX
       for ny=1:NGY
           xx(nx,ny)= a1(1) * n1(nx) + a2(1) * n2(ny);
           yy(nx,ny)= a1(2) * n1(nx) + a2(2) * n2(ny);
       end
   end
 for nx=1:NGX*repeat(1)
     for ny=1:NGY*repeat(2)
           xx(nx,ny)= a1(1) * n1(nx) + a2(1) * n2(ny);
           yy(nx,ny)= a1(2) * n1(nx) + a2(2) * n2(ny);      
           aa=nx-NGX* floor(nx/NGX);
           bb=ny-NGY* floor(ny/NGY);
           if (aa==0) aa=NGX;end
           if (bb==0) bb=NGY; end
           zz(nx,ny)=chgSlice( aa,bb);
     end
 end
%  zz=log(zz);
   surf(xx,yy,log(zz)); shading interp; view(0,90);
   axis tight; box on; colormap hot; colorbar; axis off; axis equal;
   save data xx yy zz NGX NGY NGZ chg nz

%    figure;
%    zz=log(zz);
%    zmin=min(min(zz)); zmax=max(max(zz));
%    V=linspace(zmin, zmax,500);
%    [c,h]=contourf(xx,yy,zz); set(h,'Color','none')   
%    colormap hot; colorbar; axis off; axis equal;
 
 
function chargePlot(nanoFilePath)
spin=1; % 1 or 2

% load 'f:\download\nanodcalObject.mat';
% nanoFilePath='f:\download\nanodcalObject.mat';
% file='nanodcalObject.mat';
load(nanoFilePath);

writePOSCAR=1;
[atomType, atomPos, a]=extractPOSCAR(nanoFilePath, writePOSCAR);

% charge density.  5 numbers in a line
% rho is written in sequence of z-y-x

copyfile POSCAR CHGCAR
fid=fopen('CHGCAR','a'); % append

tmp=object.RealSpaceCharge.gridNumbers;
NGX=tmp(1); NGY=tmp(2); NGZ=tmp(3);
rho=object.RealSpaceCharge.functionValues{spin};
fprintf(fid,'\n');
fprintf(fid,'%d   %d   %d\n',[NGX NGY NGZ]);

count=0;
all=0;
for ii=1:NGZ
    for jj=1:NGY
       for kk= 1:NGX
          all=all+1;
          count=count+1;
     
          if (count>5) count = 1; end
          
          if( count==5)       
              fprintf(fid,'%f   \n', rho(kk,jj,ii));          
          else              
              fprintf(fid,'%f   ', rho(kk,jj,ii));              
          end
        end
    end
end
display('Successfully written: CHGCAR.');

    