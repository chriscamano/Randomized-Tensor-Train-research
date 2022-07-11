function BTN=AwithQ(ATN,QTN)
% A*Q
% sum 3rd mode of ATN with 2nd mode of Q
BTN.n=ATN.n;
BTN.n(:,3)=QTN.n(:,3);
BTN.n(:,1)=ATN.n(:,1).*QTN.n(:,1);
BTN.n(:,4)=ATN.n(:,4).*QTN.n(:,4);

%% First core separate
tempa=reshape(ATN.core{1},[prod(ATN.n(1,2:3)),ATN.n(1,4)]);
tempa=tempa';
tempa=reshape(tempa,[ATN.n(1,4)*ATN.n(1,2),ATN.n(1,3)]);
tempq=reshape(QTN.core{1},[QTN.n(1,2),prod(QTN.n(1,3:end))]); 
B=tempa*tempq;
B=reshape(B,[ATN.n(1,4),ATN.n(1,2),QTN.n(1,3:4)]);
B=permute(B,[2,3,4,1]);
% B=reshape(B,[ATN.n(1,4),ATN.n(1,2)*prod(QTN.n(1,3:4))]);
% B=B';
BTN.core{1}=reshape(B,[ATN.n(1,2),QTN.n(1,3),QTN.n(1,4)*ATN.n(1,4)]);
for i=2:size(ATN.n,1)
    tempa=reshape(ATN.core{i},ATN.n(i,:));
    tempa=permute(tempa,[1,2,4,3]);
    tempa=reshape(tempa,[prod(ATN.n(i,1:2))*ATN.n(i,4),ATN.n(i,3)]);
    tempq=reshape(QTN.core{i},[QTN.n(i,1),prod(QTN.n(i,2:end))]); 
    tempq=tempq';
    tempq=reshape(tempq,[QTN.n(i,2),QTN.n(i,4)*QTN.n(i,1)]);    
    B=tempa*tempq;
    B=reshape(B,[ATN.n(i,[1,2,4]),QTN.n(i,[4,1])]);
    B=permute(B,[5,1,2,4,3]);
    
%     B=reshape(B,[prod(ATN.n(i,[1,2,4]))*QTN.n(i,4),QTN.n(i,1)]);
%     B=B';    
%     B=reshape(B,[QTN.n(i,1)*ATN.n(i,1)*ATN.n(i,2),ATN.n(i,4),QTN.n(i,4)]);
%     B=permute(B,[1,3,2]);
    BTN.core{i}=reshape(B,[QTN.n(i,1)*ATN.n(i,1),ATN.n(i,2),ATN.n(i,4)*QTN.n(i,4)]);
end
end