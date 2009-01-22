load ('C:/Documents and Settings/juribe/My Documents/network-project-code/sparseMatrixV1.txt');
A=spconvert(sparseMatrixV1);
[a] = textread('C:/Documents and Settings/juribe/My Documents/network-project-code/bVector.txt', '', 'delimiter', ',');
b = A\a;
save 'C:/Documents and Settings/juribe/My Documents/network-project-code/b.txt' b -ascii;