clear all
load K1.txt
load K2.txt
load K3.txt


    A=1;

    if A == 1 
        K1=log10(K1);
        K2=log10(K2);
        K3=log10(K3);
    end


%fid = fopen('three_layers_adjoint.bpp.1_6','r');
%title = fgetl(fid);
%PARS = fscanf(fid,'%*s %*s %*d %g');

for ck = 1:6

    fid = fopen(['case' num2str(ck) '.bpp.fin'],'r');
    title = fgetl(fid);
    PARS = fscanf(fid,'%*s %*s %*d %g %*g %*g');
    close all
    cfig = figure;
    if A == 1 
        PARS=log10(PARS);
    end
    PARS1=reshape(PARS(1:1400),35,40)';
    PARS2=reshape(PARS(1401:2800),35,40)';
    PARS3=reshape(PARS(2801:4200),35,40)';

    Min1=min(min(K1));
    Min2=min(min(K2));
    Min3=min(min(K3));
    MinK = min([Min1;Min2;Min3]);
    Max1=max(max(K1));
    Max2=max(max(K2));
    Max3=max(max(K3));
    MaxK = max([Max1;Max2;Max3]);



    subplot('position', [0.05 0.67 0.4 0.3])
    imagesc(K1,[MinK-1e10*eps, MaxK+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.05 0.35 0.4 0.3])
    imagesc(K2,[MinK-1e10*eps, MaxK+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.05 0.03 0.4 0.3])
    imagesc(K3,[MinK-1e10*eps, MaxK+1e10*eps]);colorbar
    axis('equal')


    subplot('position', [0.55 0.67 0.4 0.3])
    imagesc(PARS1,[MinK-1e10*eps, MaxK+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.55 0.35 0.4 0.3])
    imagesc(PARS2,[MinK-1e10*eps, MaxK+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.55 0.03 0.4 0.3])
    imagesc(PARS3,[MinK-1e10*eps, MaxK+1e10*eps]);colorbar
    axis('equal')
    print(cfig,'-depsc',['case' num2str(ck) 'pars.eps']);
    
    
    cfig2 = figure;

    subplot('position', [0.05 0.67 0.4 0.3])
    imagesc(K1,[Min1-1e10*eps, Max1+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.05 0.35 0.4 0.3])
    imagesc(K2,[Min2-1e10*eps, Max2+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.05 0.03 0.4 0.3])
    imagesc(K3,[Min3-1e10*eps, Max3+1e10*eps]);colorbar
    axis('equal')


    subplot('position', [0.55 0.67 0.4 0.3])
    imagesc(PARS1,[Min1-1e10*eps, Max1+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.55 0.35 0.4 0.3])
    imagesc(PARS2,[Min2-1e10*eps, Max2+1e10*eps]);colorbar
    axis('equal')
    subplot('position',[0.55 0.03 0.4 0.3])
    imagesc(PARS3,[Min3-1e10*eps, Max3+1e10*eps]);colorbar
    axis('equal')
    print(cfig,'-depsc',['case' num2str(ck) 'pars.eps']);
        print(cfig2,'-depsc',['case' num2str(ck) '.eps']);
    
end