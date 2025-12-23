function latexList = make_latex_examples(spaceCAS)

%
COUNT_POLY  = 0;
COUNT_RAT   = 0;
COUNT_IRRAT = 0;
%
latexList   = [];
%latexList   = '\small ';
latexList   = [latexList '\begin{center} \begin{longtable}{|c|l|l|l|} \hline '];
latexList   = [latexList 'Case & Ref. & Information & Function \\ \hline\hline '];
for CAS = spaceCAS
    COLOR       = '';
    [H,info]    = mlf.examples(CAS);
    for i = 1:numel(info.tag)
        if strcmp(info.tag{i},'polynomial')
            COUNT_POLY  = COUNT_POLY+1;
            COLOR       = '\blue';
        end
        if strcmp(info.tag{i},'rational')
            COUNT_RAT = COUNT_RAT+1;
        end
        if strcmp(info.tag{i},'irrational')
            COUNT_IRRAT = COUNT_IRRAT+1;
            COLOR       = '\orange';
        end
    end
    latexList = [latexList ...
                 COLOR '{\#' num2str(CAS) '} & ' ... 
                 info.cite ' & ' ...
                 ' $\ord= ' num2str(info.n) '$, ' info.tensor_names ' & ' ... 
                 info.name ' \\' ];
end
latexList = [latexList '\hline \caption{List of examples. \blue{Polynomial: ' num2str(COUNT_POLY) '}, {rational: ' num2str(COUNT_RAT) '}, \orange{irrational: ' num2str(COUNT_IRRAT) '}.} \label{tab:examples} \end{longtable} \end{center}\normalsize' ];
