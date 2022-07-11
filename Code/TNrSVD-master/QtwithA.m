function BTN=QtwithA(QTN,ATN)
% Q'*A
% sum 2nd mode of QTN with 2nd mode of ATN
BTN.n=ATN.n;
BTN.n(:,2)=QTN.n(:,3);
BTN.n(:,1)=ATN.n(:,1).*QTN.n(:,1);
BTN.n(:,4)=ATN.n(:,4).*QTN.n(:,4);

%% First core separate
tempq=reshape(QTN.core{1},QTN.n(1,2:4));
tempq=permute(tempq,[2,1,3]);
tempq=reshape(tempq,[prod(QTN.n(1,2:3)),QTN.n(1,4)]);
tempq=tempq';
tempq=reshape(tempq,[QTN.n(1,4)*QTN.n(1,3),QTN.n(1,2)]);
tempa=reshape(ATN.core{1},[ATN.n(1,2),prod(ATN.n(1,3:end))]); 
B=tempq*tempa;
B=reshape(B,[QTN.n(1,4),QTN.n(1,3)*prod(ATN.n(1,3:4))]);
B=B';
BTN.core{1}=reshape(B,[QTN.n(1,3),ATN.n(1,3),ATN.n(1,4)*QTN.n(1,4)]);
for i=2:size(ATN.n,1)
    tempq=reshape(QTN.core{i},[QTN.n(i,1),prod(QTN.n(i,2:end))]); 
    tempq=tempq';
    tempq=reshape(tempq,[QTN.n(i,2),QTN.n(i,4)*QTN.n(i,1)]);    
    tempa=reshape(ATN.core{i},ATN.n(i,:));
    tempa=permute(tempa,[1,3,4,2]);
    tempa=reshape(tempa,[prod(ATN.n(i,[1,3,4])),ATN.n(i,2)]);
    B=tempa*tempq;
    B=reshape(B,[(ATN.n(i,[1,3,4])),QTN.n(i,4),QTN.n(i,1)]);
    B=permute(B,[1,5,2,3,4]);
%     B=reshape(B,[prod(ATN.n(i,[1,3,4]))*QTN.n(i,4),QTN.n(i,1)]);
%     B=B';    
%     B=reshape(B,[QTN.n(i,1),ATN.n(i,1),ATN.n(i,3)*ATN.n(i,4)*QTN.n(i,4)]);
%     B=permute(B,[2,1,3]);
    BTN.core{i}=reshape(B,[QTN.n(i,1)*ATN.n(i,1),ATN.n(i,3),ATN.n(i,4)*QTN.n(i,4)]);
end

end