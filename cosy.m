j=1;
for i=1:496
    if strcmp(textdata(i,4), 'H')
        list(j,1)=textdata(i,2); list(j,2)=textdata(i,3);
        cosy(j,1)=data(i,1); cosy(j,2)=data(i+1,1);
        j=j+1;
        if strcmp(textdata(i+1,4), 'HA2')
            list(j,1)=textdata(i,2); list(j,2)=textdata(i,3);
            cosy(j,1)=data(i,1); cosy(j,2)=data(i+2,1);
            j=j+1;
        end
    end
end

scatter(cosy(:,1), cosy(:,2))