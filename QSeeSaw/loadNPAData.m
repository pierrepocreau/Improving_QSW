% Loads the hierarchy result from data files

v0plusv1 = 2;

QSW_NPA_C3 = load('data/NPA_C3.txt')';
QSW_NPA_C5_00 = load('data/NPA_C5_00.txt')';
QSW_NPA_C5_01 = load('data/NPA_C5_01.txt')';

n_C3 = length(QSW_NPA_C3);
n_C5_00 = length(QSW_NPA_C5_00);
n_C5_01 = length(QSW_NPA_C5_01);

v0_C3 = 0:2/(n_C3-1):2;
v0_C5_00 = 0:2/(n_C5_00-1):2;
v0_C5_01 = 0:2/(n_C5_01-1):2;

v1_C3 = v0plusv1*ones(1,n_C3) - v0_C3;
v1_C5_00 = v0plusv1*ones(1,n_C5_00) - v0_C5_00;
v1_C5_01 = v0plusv1*ones(1,n_C5_01) - v0_C5_01;

ratio_C3 = v0_C3./v1_C3;
ratio_C5_00 = v0_C5_00./v1_C5_00;
ratio_C5_01 = v0_C5_01./v1_C5_01;