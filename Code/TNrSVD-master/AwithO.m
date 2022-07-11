function YTN=AwithO(ATN,OTN);
% sum 3rd mode of ATN with 2nd mode of O
YTN.n=ATN.n;
YTN.n(:,3)=OTN.n(:,3);

%% First core separate
tempa=reshape(ATN.core{1},[prod(ATN.n(1,2:3)),ATN.n(1,4)]);
tempa=tempa';
tempa=reshape(tempa,[ATN.n(1,4)*ATN.n(1,2),ATN.n(1,3)]);
tempo=reshape(OTN.core{1},[OTN.n(1,2),OTN.n(1,3)]); % r(1)=r(4)=1
Y=tempa*tempo;
Y=reshape(Y,[ATN.n(1,4),ATN.n(1,2)*OTN.n(1,3)]);
Y=Y';
YTN.core{1}=reshape(Y,[ATN.n(1,2),OTN.n(1,3),ATN.n(1,4)]);
for i=2:size(ATN.n,1)
    tempa=reshape(ATN.core{i},[prod(ATN.n(i,1:3)),ATN.n(i,4)]);
    tempa=tempa';
    tempa=reshape(tempa,[ATN.n(i,4)*prod(ATN.n(i,1:2)),ATN.n(i,3)]);
    tempo=reshape(OTN.core{i},[OTN.n(i,2),OTN.n(i,3)]); % r(1)=r(4)=1
    Y=tempa*tempo;
    Y=reshape(Y,[ATN.n(i,4),ATN.n(i,1)*ATN.n(i,2)]);
    Y=Y';
    YTN.core{i}=reshape(Y,[ATN.n(i,1:2),OTN.n(i,3),ATN.n(i,4)]);
end
end