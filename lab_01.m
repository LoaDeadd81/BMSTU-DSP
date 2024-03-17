sigma = 5;
plot_step = 0.005;

n = input('Input n: ');
dt = input('Input dt: ');

t_max = dt*(n-1)/2;
t = -t_max:dt:t_max; 
l = -2;
r = 2;

rect_discrete = rect_arr(t, l, r);
gauss_discrete = gauss_arr(t, sigma);

x = -t_max:plot_step:t_max;
gauss_restored = zeros(1, length(x));
rect_restored = zeros(1, length(x));
for i=1:length(x)
   for j = 1:n
       tmp = (x(i)-t(j)) * pi / dt;
       gauss_restored(i) = gauss_restored(i) + gauss_discrete(j) * sinc(tmp);
       rect_restored(i) = rect_restored(i) + rect_discrete(j) * sinc(tmp);
   end
end

rect_ref = rect_arr(x, l, r);
gauss_ref = gauss_arr(x, sigma);

figure;

subplot(2,1,1);
title('Прямоугольный импульс');
hold on;
grid on;
plot(x, rect_ref, 'k');
plot(x, rect_restored, 'b');
plot(t, rect_discrete, '.g');
legend('Исходная', 'Восстановленная', 'Дискретная');

subplot(2,1,2);
title('Гауссовский фильтр');
hold on;
grid on;
plot(x, gauss_ref, 'k');
plot(x, gauss_restored, 'b');
plot(t, gauss_discrete, '.m');
legend('Исходная', 'Восстановленная', 'Дискретная');

function y = sinc(x)
    y = sin(x)/x;
end

function arr = rect_arr(x, l, r)
    arr = zeros(size(x));
    for i=1:length(x)
        arr(i) = rect_f(x(i), l, r);
    end
end

function y = rect_f(x, l, r)
    if l <= x && x <= r
        y = 1;
    else
        y = 0;
    end
end

function arr = gauss_arr(x, sigma)
    arr = zeros(size(x));
    for i=1:length(x)
        arr(i) = gauss_f(x(i), sigma);
    end
end

function y = gauss_f(x, sigma)
    y = exp(-x^2/sigma^2);
end
