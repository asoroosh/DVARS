clear

addpath(genpath('~/.../DVARS'))
%addpath(genpath('~/.../xDF')) % Just in case you wanted to do FDR instead
%of Benjamini-Hochberg

clear
I = 1000;
T = 500;

% pure, spike less noise:
Yp = randn(I,T);

% put a corruption on datapoint 250:
Yc = [Yp(:,1:250) Yp(:,251:end)+10];

[I,T] = size(Yp); 

% Run octave version % should find a corruption on datapoint 250!
[DVARS_oct,Stat_oct]=DVARSCalc_octave(Yc);
BH_pval_adj_oct = Stat_oct.pvals.*T;
find(BH_pval_adj_oct<0.05)

% Run the matlab version % should find a corruption on datapoint 250!
[DVARS_mat,Stat_mat]=DVARSCalc(Yc);
BH_pval_adj_mat = Stat_mat.pvals.*T;
find(BH_pval_adj_mat<0.05)

% compare the p-values
figure; 
hold on; 
scatter(BH_pval_adj_mat,BH_pval_adj_oct)
ylabel('Matlab','Interpreter','latex')
xlabel('Octave','Interpreter','latex')