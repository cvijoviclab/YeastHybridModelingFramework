function [rapamycin,SNF1deletion]=loadDEletionProteomics()
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/linoste/Documents/GitHub/Crabtree/proteomics/rapamycin.xls
%    Worksheet: A. I(g) Intensity
%
% Auto-generated by MATLAB on 01-Jun-2020 19:14:00

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 17);

% Specify sheet and range
opts.Sheet = "A. I(g) Intensity";
opts.DataRange = "A1:Q4149";

% Specify column names and types
opts.VariableNames = ["ORF", "Genename", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "RAP700"];
opts.SelectedVariableNames = ["ORF", "Genename", "RAP700"];
opts.VariableTypes = ["string", "string", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "char", "double"];

% Specify variable properties
opts = setvaropts(opts, ["ORF", "Genename", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ORF", "Genename", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16"], "EmptyFieldRule", "auto");

% Import the data
rapamycin = readtable("/Users/linoste/Documents/GitHub/Crabtree/proteomics/rapamycin.xls", opts, "UseExcel", false);


%% Clear temporary variables
clear opts
%% Import data from spreadsheet
% Script for importing data from the following spreadsheet:
%
%    Workbook: /Users/linoste/Documents/GitHub/Crabtree/proteomics/SNF1deletion.xlsx
%    Worksheet: Sheet1
%
% Auto-generated by MATLAB on 01-Jun-2020 20:23:19

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 7);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A1:G365";

% Specify column names and types
opts.VariableNames = ["YGL205W", "POX1", "Var3", "Var4", "Var5", "Var6", "VarName7"];
opts.SelectedVariableNames = ["YGL205W", "POX1", "VarName7"];
opts.VariableTypes = ["string", "string", "char", "char", "char", "char", "double"];

% Specify variable properties
opts = setvaropts(opts, ["YGL205W", "POX1", "Var3", "Var4", "Var5", "Var6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["YGL205W", "POX1", "Var3", "Var4", "Var5", "Var6"], "EmptyFieldRule", "auto");

% Import the data
SNF1deletion = readtable("/Users/linoste/Documents/GitHub/Crabtree/proteomics/SNF1deletion.xlsx", opts, "UseExcel", false);


%% Clear temporary variables
clear opts