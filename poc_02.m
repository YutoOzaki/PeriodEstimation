%%
T = readtable('Consumo_cerveja.csv');
x = table2array(T(:, end));

%%
lmd = [1, 1e-3, 1e-2];
Pmax = 50;
P = PIE(x, Pmax, lmd);

%%
figure(1);
clf; cla;
subplot(2, 1, 1);
plot(x);
subplot(2, 1, 2);
stem(P);
xlim([2, Inf]);