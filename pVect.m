load sparseMatrixV1.txt;
A = spconvert(sparseMatrixV1);
[a] = textread('bVector.txt', '', 'delimiter', ',');
b = A\a;
save b.txt b -ascii;