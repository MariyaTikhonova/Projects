
# coding: utf-8

# # Алгоритм проверки того, что базис действительно является базисом Гребнера

# В Sage не было классов для реализации многочленов с переменными индексами, поэтому пришлось реализовать класс мономов с переменными индексами и полиномов. Мне привычнее работать с питоном и для реализации алгоритма я пока не использовала специфических свойств Sage, поэтому реализовала данные классы на Python (в питоновских тетрадках с ними работать удобнее). Однако данный код без труда можно перенестии в Sage, для расширения дальнейших возможностей классов.

# In[1]:

import numpy as np
import fractions


# # Реализация классов мономов и многочленов с переменными индексами.

# Я реализовала базовый класс мономов с переменными индексами. Для них я переопределила базовые операции (сложение, умножение, вывод на печать, сравнение (у меня оно идет в лексикографическом порядке), проверку на равенство). Также добавила возможность замены переменного индекса с помощью функции m.change_var('k', 'l'). Для мономов определила операцию сдвига индексов с помощью операции m.make_shift('k', -1), которая сдвигает соответственной все индексы по переменной на заданное число позиций и m.normalize_inds(), которая делает сдвиги по каждому переменному индексу, чтобы минимальный (или максимальный из них имел сдвиг 0). 

# Также я реализовала деление одного монома на другой, m.calc_devisor(delimiter). Правда пока с учетом названий переменных индексов. Результатом будет моном с нормализованными индексами.

# Для классов полиномов переопределила базывые функции, функции pol.lt(), аналогичные функции сдвигов и нормализации (нормализация производится относительно lt). Также реализовала функции проверки на однородность по степени и на то, что многочлен нулевой.

# Ниже будут приведены примеры задания и использования данных классов.

# In[2]:


class Monom:
    def __init__(self, fixed_part = {}, var_part = {}, koef = 1):
        self.fixed_part = {}
        self.var_part = {}
        self.koef = koef
        for fixed_ind in fixed_part.keys():
            if fixed_part[fixed_ind] != 0:
                self.fixed_part[fixed_ind] = fixed_part[fixed_ind]
        for var in var_part.keys():
            self.var_part[var] = {}
            for var_ind in var_part[var].keys():
                if var_part[var][var_ind] != 0:
                    self.var_part[var][var_ind] = var_part[var][var_ind]
            if len(self.var_part[var].keys()) == 0:
                self.var_part.pop(var)

    def __neg__(self):
        return -1 * self
    
    def __add__(self, other):
        if self == other:
            return Monom(self.fixed_part, self.var_part, self.koef + other.koef) 
        return Monom()
    
    def __mul__(self, other):
        ##print type(other)
        if type(other) in (type(0.5), type(1), type(-1), type(-0.1)):
            return Monom(self.fixed_part, self.var_part, other * self.koef)
        
        new_fixed_part = {}
        new_var_part = {}
        for f in self.fixed_part:
            new_fixed_part[f] = self.fixed_part[f]
        for f2 in other.fixed_part:
            if f2 in new_fixed_part.keys():
                new_fixed_part[f2] += other.fixed_part[f2]
            else: 
                new_fixed_part[f2] = other.fixed_part[f2]
        
        for var in self.var_part:
            new_var_part[var] = {}
            for var_ind in self.var_part[var]:
                new_var_part[var][var_ind] = self.var_part[var][var_ind]
        for var2 in other.var_part:
            if var2 not in new_var_part.keys():
                new_var_part[var2] = {}
            for var_ind2 in other.var_part[var2]:
                if var_ind2 not in new_var_part[var2].keys():
                    new_var_part[var2][var_ind2] = other.var_part[var2][var_ind2]
                else:
                    new_var_part[var2][var_ind2] += other.var_part[var2][var_ind2]                
        return Monom(new_fixed_part, new_var_part, self.koef * other.koef)
    
    def __rmul__(self, other): 
        return self.__mul__(other)
    
    def __sub__(self, other):
        p = -1 * other
        return self + p
    
    def __lt__(self, other):
        flag = False
        if self.fixed_part < other.fixed_part:
            flag = True
        if len(self.var_part) == 0 and len(other.var_part) == 0:
            return flag
        if len(self.var_part) == 0 and len(other.var_part) != 0:
            return True
        if len(self.var_part) != 0 and len(other.var_part) == 0:
            return False
        return self.var_part < other.var_part
    

    def  __eq__(self, other):
        return self.fixed_part == other.fixed_part and self.var_part == other.var_part
    
    def  __gt__(self, other):
        return not (self < other) and not (self == other) 
    
    def  __ge__(self, other):
        return not (self < other)
            
    def __str__(self):
        string = ''
        if abs(self.koef) < 0.00001:
            return str(0)
        if abs(self.koef - 1) > 0.0001:
            string += str(self.koef) + ' '
            
        ##print(self.fixed_part.keys())
        for fixed_ind in self.fixed_part.keys():
            if self.fixed_part[fixed_ind] > 0:
                string += 'x_'+str(fixed_ind)
                if self.fixed_part[fixed_ind] > 1:
                    string += '^('+str(self.fixed_part[fixed_ind])+')'
                ##print('x_'+str(fixed_ind)+'^('+str(fixed_part[fixed_ind])+')')
                string += ' '
        for var in self.var_part.keys():
            for var_ind in self.var_part[var].keys():
                if self.var_part[var][var_ind] > 0:
                    string += 'x_('+ var
                    if var_ind > 0 :
                        string += '+'+ str(var_ind)
                    if var_ind < 0:
                        string += str(var_ind)
                    string += ')'
                if self.var_part[var][var_ind] > 1:
                    string += '^('+str(self.var_part[var][var_ind])+')'
                string += ' '
        if len(string) == 0:
            return str(self.koef)
        return string
    
    def __repr__(self):
        return self.__str__()

    def change_var(self, a, b):
        if a in self.var_part.keys() and not b in self.var_part.keys():
            tmp = self.var_part.pop(a)
            self.var_part[b] = tmp
            
    def make_shift(self, key, shift):
        new_dict = {}
        if key in self.var_part.keys():
            for ind in self.var_part[key]:
                new_dict[ind + shift] = self.var_part[key][ind]
            self.var_part[key] = new_dict.copy()
    
    def normalize_shift(self, key, rev = False):
        max_shift = max(self.var_part[key])
        if rev:
            max_shift = min(self.var_part[key])
        self.make_shift(key, -max_shift)
    
    def normalize_inds(self, rev = False):
        for key in self.var_part:
            self.normalize_shift(key, rev)
## Помни, что здесь вызавается нормализация => меняются исходные мономы            
    def calc_devisor(self, other, change = False, norm = True):
        if not change:
            m1 = Monom(self.fixed_part, self.var_part)
            m2 = Monom(other.fixed_part, other.var_part)
            return m1.calc_devisor(m2, change = True, norm = norm)
        new_fixed_part = {}
        new_var_part = {}
        if norm:
            self.normalize_inds()
            other.normalize_inds()
        #print(self, other)
        for ind in other.fixed_part:
            if not self.fixed_part.get(ind):
                return None
            power = self.fixed_part[ind] - other.fixed_part[ind]
            if power < 0:
                return None
            if power > 0:
                new_fixed_part[ind] = power
        for ind in self.fixed_part:
            if not other.fixed_part.get(ind):
                new_fixed_part[ind] = self.fixed_part[ind]
                
        #print(self, other)
        for var in other.var_part: 
            if not self.var_part.get(var):
                return None
            new_var_part[var] = {}
            for ind in other.var_part[var]:
                if not self.var_part[var].get(ind):
                    return None
                power = self.var_part[var][ind] - other.var_part[var][ind]
                if power < 0:
                    return None
                if power > 0:
                    new_var_part[var][ind] = power
            for ind in self.var_part[var]:
                if not other.var_part[var].get(ind):
                    new_var_part[var][ind] = self.var_part[var][ind]
            if len(new_var_part[var].keys()) == 0:
                new_var_part.pop(var)
        #print(new_fixed_part)
        for ind in self.var_part:
                if not other.var_part.get(ind):
                    new_var_part[ind] = self.var_part[ind].copy()
        m = Monom(new_fixed_part, new_var_part)
        return m
    
    def calc_power(self):
        power = np.sum(self.fixed_part.values())
        for var in self.var_part:
            power += np.sum(self.var_part[var].values())
        return power
            


            
class Polynom:
    def __init__(self, monoms):
        self.monoms = []
        for mm in monoms:
            m = Monom(mm.fixed_part, mm.var_part, mm.koef)
            if m in self.monoms:
                ind = self.monoms.index(m)
                self.monoms[ind] += m
            else:
                self.monoms.append(m)
        self.check()
                
    def __str__(self):
        string = ''
        flag = False
        for monom in self.monoms:
            if flag and monom.koef >= 0:
                string += '+ '
            flag = True
            string += monom.__str__()
        if len(string) == 0:
            return '0'
        return string
    
    def __repr__(self):
        return self.__str__()
    
    def __add__(self, other):
        new_monoms = self.monoms + other.monoms
        self.check()
        return Polynom(new_monoms)
    
    def __iadd__(self, other):
        p = self + other
        self.monoms = p.monoms
        self.check()
        return self
    
    def __sub__(self, other):
        p = -1 * other
        self.check()
        return self + p
    
    def __isub__(self, other):
        p = self - other
        self.monoms = p.monoms
        return self
    
    def __mul__(self, other):
        new_monoms = []
        ##print '!!', other
        ##if type(other) in (float, int):
        for  m in self.monoms:
            new_monoms.append(other * m)
        return Polynom(new_monoms)
    
    def __rmul__(self, other): 
        return self.__mul__(other)
    
    def __neg__(self):
        return -1 * self
    
    def lt(self):
        return max(self.monoms)
    
    def change_var(self, a, b):
        for m in self.monoms:
            m.change_var(a,b)
            
    def make_shift(self, key, shift):
        for m in self.monoms:
            m.make_shift(key, shift)
    
    def normalize(self, rev = False):
        s_l = self.lt()
        for var in s_l.var_part.keys():
            if rev:
                shift = min(s_l.var_part[var].keys())
            else:
                shift = max(s_l.var_part[var].keys())
            for m in self.monoms:
                m.make_shift(var, -shift)
        
    
    def check(self):
        for m in self.monoms:
            if abs(m.koef) < 0.00001:
                #a = 6
                self.monoms.pop(self.monoms.index(m))
    def copy(self):
        return Polynom(self.monoms)
    
    def is_zero(self):
        self.check()
        if len(self.monoms) == 0:
            return True
        else:
            return False
    
    def is_uniform(self):
        if len(self.monoms) == 0:
            return True
        power = self.monoms[0].calc_power()
        for m in self.monoms:
            if not m.calc_power() == power:
                return False
        return True


# # Примеры использования данных классов.

# In[3]:

#Моном задается с помощью словарей.
#Словарь фиксированных индексов: key - индекс, value - степень
fixed_part = {3 : 1, 4 : 5}
#Переменные индексы задаются с помощью словаря словарей:
#Индекс переменной, а для нее определяется следующий словарь:
# key - сдвиг, value - степень
var_part = {'k' : {0 : 1, -2 : 2}, 'l' : {0 : 3}}
# Коэфиициент при мономе. По умолчанию 1.
koef = 6
m = Monom(fixed_part, var_part, koef)
m


# In[4]:

m.change_var("l", "s")
print m, '\n'
m.make_shift('k', 100)
print m, '\n'
m.normalize_inds()
print m, '\n'


# Перед делением мы нормализуем индексы, поэтому получается следующий результат.

# In[5]:

n = Monom({4 : 2}, {'k' : {100 : 1}})
print 'm = ', m , '\n n = ', n, '\n m/n = ', m.calc_devisor(n)


# Можно убрать нормализацию.

# In[6]:

print 'm = ', m , '\n n = ', n, '\n m/n = ', m.calc_devisor(n, norm = False)


# При разных переменных индексах однако деление не выполняется.

# In[7]:

n.change_var('k', 't')
print 'm = ', m , '\n n = ', n, '\n m/n = ', m.calc_devisor(n)


# In[8]:

#Задаем полином в виде списка мономов.
pol = Polynom([m, n])
pol


# Сдвиги и нормализация (производится по старшему члену).

# In[9]:

pol.make_shift("k", 5)
pol.make_shift("t", -95)
pol


# In[10]:

pol.normalize()
pol


# In[11]:

pol2 = Polynom([-5*m])
pol2


# Полиномы можно складывать, умножать на константы и мономы (правда умножать на мономы нужно справа, у меня возникли трудности при переопределении оператора для умножения моном * полином)

# In[12]:

print - 7 * pol - 1/2. * pol2
print - 1/30. * pol2 * m 


# Проверка, что многочлен нулевой. А также проверка того, что он однороден по степени.

# In[13]:

print pol.is_zero()
print (pol - pol).is_zero()
print pol.is_uniform()


# In[14]:

#Копирование лучше осуществлять с помощью функции copy().
f = pol.copy()
f


# # Алгоритм проверки

# Алгоритм проверки, является ли данный базис базисом Гребнера своего идеала. Алгоритм работает для различных базисов, содержащих более одного многочлена, или для неавторедуцированных базисов. Однако в данных случаях он может зацикливаться, поэтому я установила ограничение на число операций при редукции.
# Например, он, как мы с Вами обсуждали, зацикливается на следующем базисе:
# g = x_(k) - x_(k-1). Далее я произвожу тестирование алгоритма. Это будет видно на примерах.

# In[15]:

def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)


# In[16]:

def calc_gcd (m1, m2):
    new_fixed_part = {}
    new_var_part = {}
    for key in m1.fixed_part.keys() + m2.fixed_part.keys():
        if key in m1.fixed_part.keys() and key in m2.fixed_part.keys():
            new_fixed_part[key] = min(m1.fixed_part[key], m2.fixed_part[key])
    for var in m2.var_part.keys() + m2.var_part.keys():
        if var in m1.var_part.keys() and var in m2.var_part.keys():
            new_var_part[var] = {}
            for key in m1.var_part[var].keys() and m2.var_part[var].keys():
                new_var_part[var][key] = min(m1.var_part[var][key], m2.var_part[var][key])
    return Monom(new_fixed_part, new_var_part, fractions.gcd(m1.koef, m2.koef))
                
def calc_S_pair(f, g):
    gcd = calc_gcd(f.lt(), g.lt())
    m = (f.lt() * g.lt()).calc_devisor(gcd) 

    m1 = m.calc_devisor(g.lt())
    m2 = m.calc_devisor(f.lt())
    #print m1, m2, f.lt(), g.lt()
    S = f * m2 - g * m1
    #in some cases the sign was not correct
    if not S.is_zero() and S.lt() == m2 * f.lt():
        S_2 = f * m2 + g * m1 
        if S_2.lt() != (f * m2).lt():
            #print S, S_2, m2 * f.lt(), m1 * g.lt(), f * m2, g * m1 
            return S_2
    return f * m2 - g * m1



# In[17]:

def autoreduce_check(G):
    lts = []
    for g in G:
        lts.append(g.lt())
    for i in range(len(lts)):
        lt = lts[i]
        for g in G:
            for m in g.monoms:
                if not m == lt:
                    if m.calc_devisor(lt):
                    #print m, lt
                        return False
        for j in range(len(lts)):
            if not i == j:
                if lts[j].calc_devisor(lt):
                    return False
    return True
            


# In[18]:

def reduce_S_pol(S, G, iter = 0, shift = 0):
    for k in range(len(G)):
        g = G[(k + iter + shift) % len(G)].copy()
        ##print "!!!!", G
        #.normalize()
        # print 'S = ', S, "##", g #,  len(G), iter, (k + iter + shift) % len(G)
        g.normalize()
        s_l = S.lt()
        g_l = g.lt()
    
        m = s_l.calc_devisor(g_l)
        if not m:
            ##print "(((("
            continue
        # print "Reduction  ",  'S = ', S, "##", g, "S_lt =  ", S.lt()
        red =  g * m
    #print s_l, m, red
        for var in s_l.var_part:
            shift = min(s_l.var_part[var])
            red.make_shift(var, shift)
        
    #print s_l, m, red
        return red.lt().koef * S - s_l.koef * red, True
    return S, False

def reduce_check_S_pol(S, G, show = False, iter = 5):
    show = False
    if S.is_zero():
        return True
    #print S, S.lt(), S.is_zero(), 
    i = 0
    S, flag = reduce_S_pol(S, G, i)
    
    iter = 5
    if show:
        print S, S.lt(), S.is_zero(), 
    while flag and i < iter and not S.is_zero():
        S, flag = reduce_S_pol(S, G, i)
        #S.normalize()
        if show:
            print S
        i += 1
    #if i == iter:
        #print "Iteration limit exceeded!"
    if S.is_zero():
        return True
    return False

def basis_check(G, show = False, max_comb = 3):
    max_len = 1
    G_new = []
    if len(G) > max_len:
        print "Too many Polynomials in the basis, the algorithm may not stop.\n"
        #return False
    S_set = set()
    for g in G:
        for g_ in G:
            g1 = g.copy()
            g2 = g_.copy()
            g1.normalize()
            g2.normalize()
            S_set.add(calc_S_pair(g1, g2))
            G_new.append(g1)
            G_new.append(g2)
            for i in range(len(g2.lt().var_part.keys())):
                for keys in combinations(g2.lt().var_part.keys(), i + 1):
        #print keys
                    for k in keys:
            #print key
            #for k in g2.lt().var_part.keys():
                        if k in g1.lt().var_part.keys():
                            g2.change_var(k, k + k)
                            S_set.add(calc_S_pair(g1, g2))
                            G_new.append(g2)
    G_new = list(set(G_new))
    while len(S_set) > 0:
        #print "Reduction try!"
        S = S_set.pop()
        if not reduce_check_S_pol(S.copy(), G_new, show):
            
            print "Extented algorithm is in process this may take some time.\n"
            if not reduce_dif_S_pol(S.copy(), G_new, max_comb):
                print "This S_pair could not be reduced: ", S #, "    ", S.lt(), '\n'
            #S_save = S
                return S, False
            print "Extended algorithm reduced this S-pair: ", S, "\n"
    return Polynom([]), len(S_set) == 0


# После тестирование я определила, что зацикливание происходит из-за того, что мы пытаемся отредуцировать многочлены по базису в неправильном порядке (в алгоритме результат редукции сущестенно зависит от порядка многочленов в базисе). Поэтому в случае если у нас многочлен не удается отредуцировать по базису стандартным способом ('iteration limit exceeded'), теперь запускаестся расширенный алгоритм, который пробует отредуцировать многочлен по различным перестановкам базиса.

# Его идея состоит в следующем:  пройтись по всевозможным перестановкам базиса и для каждой попытаться отредуцировать многочлен по данной перестановке. 

# Чтобы избежать вычислительных трудностей (количество перестановок растет как факториал от размера базиса, а по построению в базисе у меня находятся всевозможные многочлены со всевозможными значениями индексов), я ввела параметр max_comb. Вместо всевозможных перестановок мы проходимся по всевозможным подмножествам элементов базиса длины max_comb и для каждого из них пытаемся отредуцировать многочлен по всевозможным перестановкам данного подмножества.
# 

# По умолчанию max_comb равен 3.

# In[19]:

import itertools
def reduce_dif_S_pol(S_old, G, r = 3):
    for G_comb in combinations(G, r):
        #if reduce_check_S_pol(S.copy(), G_comb): 
        #    return True
        for G_perm in itertools.permutations(G_comb):
            S = S_old.copy()
            #print "!!"
            #print G_perm
            for i in range(len(G_perm)):
                S, flag = reduce_S_pol(S, [G_perm[i]])
                #print S
                if S.is_zero():
                    return True
        
                    
            

    print "All permutations have been checked!"
    return False
        


# In[20]:

def check_G(G, show = False, max_comb = 3, generation = False):
    print "Polynomials in basis: ", len(G), '\n'
    for g in G:
        print g
    print '\n'
    if not autoreduce_check(G):
        print "Basis is not autoreduced. The algorithm may not stop.\n"
    for g in G:
        if not g.is_uniform():
            print "Basis is not uniform.\n"
            break
    S, flag = basis_check(G, show, max_comb)
    print "Verification\n ", flag, "\n"
    return S


# # Тестирование алгоритма.

# Тестирования алгоритма. На вход подается множество многочленов с переменными индексами. Фукция также проверяет выполнение различных условий для базиса (авторедуцированность и т. д.). В случае если программа заканчивает работу в связи с превышением лимита операций, печатается результат False, а также сообщение о том, что лимит операций превышен.

# Ниже приведены различные тесты работы чекера.

# In[21]:

mm1 = Monom({1 : 1}, {'k': {0 : 1}})
mm2 = Monom({}, {'k': {0 : 1}})
f5 = Polynom([mm2, mm1])
check_G([f5])


# In[22]:

mm1 = Monom({1 : 1, 4 : 6, 9 : 11}, {'k': {0 : 1}})
mm2 = Monom({5 : 1}, {'k': {0 : 1}})
f5 = Polynom([mm2, mm1])
check_G([f5])


# Пример действительно зацикливается.

# In[23]:

mm1 = Monom({}, {'k': {0 : 1}})
mm2 = Monom({}, {'k': {-1 : 1}})
f5 = Polynom([mm1, mm2])
check_G([f5], show = True)


# In[24]:

mm1 = Monom({1 : 1, 4 : 6, 9 : 11}, {})
mm2 = Monom({6 : 5}, {})
f5 = Polynom([mm2, mm1])
check_G([f5])


# In[25]:

mm1 = Monom({1 : 1, 4 : 6, 9 : 11}, {})
mm2 = Monom({}, {'k': {0 : 1}})
f5 = Polynom([mm2, mm1])
check_G([f5])


# In[26]:

mm1 = Monom({1 : 1}, {'k': {0 : 1}})
mm2 = Monom({0 : 1}, {'k': {-1 : 1}})
f5 = Polynom([mm2, mm1])
check_G([f5])


# In[27]:

mm1 = Monom({1 : 1}, {'k': {0 : 1}})
mm2 = Monom({0 : 1}, {'k': {-1 : 1}})
mm3 = Monom({8 : 1}, {})
f5 = Polynom([mm2, mm1])
f6 = Polynom([mm3])
check_G([f5, f6])


# In[28]:

mm1 = Monom({1 : 1}, {'k': {0 : 2}})
mm2 = Monom({0 : 1}, {'k': {-1 : 1}})
mm3 = Monom({8 : 1}, {})
mm4 = Monom({0 : 1}, {'k': {0 : 1}})
f5 = Polynom([mm2, mm1])
f6 = Polynom([mm3, mm4])
check_G([f5, f6], max_comb = 4)


# In[29]:

mm1 = Monom({1 : 1}, {'k': {0 : 2}, 'l': {0 : 20}})
mm2 = Monom({1 : 1}, {'k': {0 : 2}})
mm3 = Monom({8 : 1}, {})
mm4 = Monom({0 : 1}, {'k': {0 : 1}})
f5 = Polynom([mm2, mm1])
check_G([f5])


# Надо быть внимательными к порядку. У нас x_1 < x_2 < x_3. Значит при проверке того, что базис является базисом Гребнера надо правильно называть переменные. Так (x + z, y - z) надо представлять так x_3 + x_1, x_2 - x_1.

# In[30]:

mm1 = Monom({1 : 1})
mm2 = Monom({2 : 1})
mm3 = Monom({3: 1})
f5 = Polynom([mm1, mm3])
f6 = Polynom([mm2, -mm1])
check_G([f5, f6])


# In[31]:

mm1 = Monom({},{'a' : {0 : 1}}) #z
mm2 = Monom({},{'k' : {0 : 1}}) #y
mm3 = Monom({},{'z' : {0: 1}}) #x
f5 = Polynom([mm1, mm3])
f6 = Polynom([mm2, -mm1])
check_G([f5, f6])


# In[32]:

# x > y > z
#Not Groenber: y - x^2, z - x^3
mm1 = Monom({1 : 1}) #z
mm2 = Monom({2 : 1}) #y
mm3 = Monom({3: 3}) #x
mm4 = Monom({3: 2}) #x
f5 = Polynom([mm1, -mm3])
f6 = Polynom([mm2, -mm4])
check_G([f5, f6])


# In[33]:

# y > z > x
#Groenber: y - x^2, z - x^3
mm1 = Monom({1 : 1}) #z
mm2 = Monom({2 : 1}) #y
mm3 = Monom({0: 3}) #x
mm4 = Monom({0: 2}) #x
f5 = Polynom([mm1, -mm3])
f6 = Polynom([mm2, -mm4])
check_G([f5, f6])


# In[34]:

mm1 = Monom({}, {'k': {0: 1}, 'l' : {1 : 1}})
mm2 = Monom({}, {'l': {0: 1}, 'k' : {1 : 1}})
f5 = Polynom([mm1, mm2])
check_G([f5], show = True)


# Однородности по весу тоже недостаточно.

# In[35]:

mm1 = Monom({}, {'k': {0 : 1}})
mm2 = Monom({1:1}, {'k': {-1 : 1}})
f5 = Polynom([mm1, mm2])
check_G([f5], show = True)


# In[36]:

mm1 = Monom({1 : 1}, {'k': {0 : 1}})
mm2 = Monom({0 : 1}, {'k': {1 : 1}})
f5 = Polynom([mm2, mm1])
check_G([f5])


# # УРА! Пример поддался!

# Для того, чтобы отредуцировать все S-полиномы в данном примере пришлось создать функцию reduce_dif_S_pol(), которая запускается, если не удается отредуцировать полином стандартным способом. Подрбный комментарий по ее работе написан над реализацией функции.

# При параментре max_comb = 1 отредуцировать все S-полиномы не удается. При max_comb >= 2 - редукция проходит успешно, однако при max_comb = 6 и больше время вычисления становится уже заметным. Для 10 я даже дожидаться не стала. Но у меня медленный компьютер.

# In[39]:

mm1 = Monom({1 : 1}, {'k': {0 : 1}})
mm2 = Monom({0 : 1}, {'k': {1 : 1}})
f5 = Polynom([mm2, mm1])
mm3 = Monom({1: 1}, {'kk': {0: 1}, 'k': {-1: 1}})
mm4 = Monom({1: 1}, {'kk': {-1: 1}, 'k': {0: 1}}, -1)
f6 = Polynom([mm3, mm4])
check_G([f5, f6], max_comb = 2)


# In[55]:

mm1 = Monom({1 : 1, 2 : 3}, {'k': {0 : 1}})
mm2 = Monom({1 : 2, 4 : 6}, {'k': {1 : 1}})
f5 = Polynom([mm2, -mm1])
mm3 = Monom({1 : 1, 2 : 3}, {'kk': {0: 1}, 'k': {-1: 1}})
mm4 = Monom({1 : 1, 2 : 3}, {'kk': {-1: 1}, 'k': {0: 1}}, -1)
f6 = Polynom([mm3, mm4])
check_G([f5, f6], max_comb = 2)


# # Есть еще куда расти!

# В самом алгориме пока остались технически нереализованные вопросы. Например, я пока еще не вводила ограничения на нижние индексы для мономов k >= k0. Это, безусловно следует сделать. 

# Я уже упоминала, что данная реализация должна работать и в Sage. Возможно надо написать оболочку, которая бы могла работать с обычными полиномами из Sage. Однако, мне кажется, разумно это сделать уже после отладки работы алгорима.

# Также буду вам благодарна, если другие недостатки в программе, я могла забыть про какие-нибудь вырожденные случаи или может подскажите, какие еще методы разумно реализовывать для классов мономов и полиномов.

# Буду Вам очень признательна за Ваши комментарии. Я долго работаю над программой, хочется довести ее до хорошего состояния.

# # Алгоритм построения базиса Гребнера

# Следующем шагом является модификация алгоритма проверки, чтобы он не просто проверял, что базис является базисом Гребнера, а по исходному множеству параметрических многочленов дополнял их до Базиса Гребнера. 

# In[56]:

def compute_basis(B, max_comb = 3, max_size = 6):
    if len(B) > max_size:
        print "Too many Polynomials in input."
        return B
    while len(B) < max_size + 1:
        s = check_G(B, max_comb)
        if s.is_zero():
            print "The basis is computed:\n"
            print B
            return B
        B.append(s)
    print "Iteration limit exceeded. The basis is computed only partially:\n"
    print B
    return B


# Алгоритм справился с модельным примером.

# In[57]:

mm1 = Monom({1 : 1}, {'k': {0 : 1}})
mm2 = Monom({0 : 1}, {'k': {1 : 1}})
f5 = Polynom([mm2, mm1])
compute_basis([f5])


# Однако с данным алгоритмом возникают огромные вычислительные трудности в связи с огромным количеством вариантов, которые приходится перебирать. Поэтому эта часть еще нуждается в дополнительной оптимизации. 

# In[ ]:



