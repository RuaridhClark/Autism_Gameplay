figure;
axis([0 1100 0 800])

for i = 1 : 16
    hold on
    scatter(obj(i,1),obj(i,3),'+')
    hold on
    scatter(obj(i,2),obj(i,4),'+')
    hold on
    scatter(obj(i,1),obj(i,4),'+')
    hold on
    scatter(obj(i,2),obj(i,3),'+')
    plot([obj(i,1),obj(i,2)],[obj(i,3),obj(i,3)],'b')
    plot([obj(i,1),obj(i,2)],[obj(i,4),obj(i,4)],'b')
    plot([obj(i,1),obj(i,1)],[obj(i,3),obj(i,4)],'b')
    plot([obj(i,2),obj(i,2)],[obj(i,3),obj(i,4)],'b')
end