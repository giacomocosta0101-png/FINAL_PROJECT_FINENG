%% Comb-Bernoulli case:

Z = randn(200,2);
mean1 = 2.56;
std1 = 1.00;
mean2 = 2.10;
std2 = 1.00;

X1 = Z(:,1);
X2 = 0.4*X1+sqrt(1-0.4)*Z(:,2);

% Combine X1 and X2 into a single matrix for further analysis
U1 = normcdf(X1);
U2 = normcdf(X2);

survival_1 = find(U1<=0.8);
survival_2 = find(U2<=0.8);

X1(survival_1) = 0;
X1(~survival_1) = exp((std1*U1(~survival_1)-0.8)./0.2+mean1);

X2(survival_2) = 0;
X2(~survival_2) = exp((std2*U2(~survival_2)-0.8)./0.2+mean2);

scatter(U1,U2)
hold on
rectangle('Position',[0 0 1 1])
axis square
axis([0 2 0 2]);
%% Semiparametric case

Z = randn(200,2);
mean1 = 2.56;
std1 = 1.00;
mean2 = 2.10;
std2 = 1.00;

X1 = Z(:,1);
X2 = 0.4*X1+sqrt(1-0.4)*Z(:,2);

% Combine X1 and X2 into a single matrix for further analysis
U1 = normcdf(X1);
U2 = normcdf(X2);

survival_1 = find(U1<=0.8);
survival_2 = find(U2<=0.8);

U1(survival_1) = 0;
U2(survival_2) = 0;

scatter(X1,U1)

%or:
%X1(survival_1) = 0;
%scatter(X1,U1)

U1 = sort(U1);
U2 = sort(U2);







