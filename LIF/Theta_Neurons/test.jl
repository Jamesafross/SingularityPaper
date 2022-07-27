mutable struct test1
    a
end

mutable struct test2
    b
end

ab = test2(ones(2))

tb = test1(b)

tb.a.b = tb.a.b*2

