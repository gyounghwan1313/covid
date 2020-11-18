import matplotlib.pylab as plt # 그래프
from matplotlib import font_manager, rc # 그래프 설정 변경
import missingno as msno # 결측치 시각화
import numpy as np
import pandas as pd
import scipy.stats as ss #z-score
import pingouin as pg # correation test
import seaborn as sns # correation table
import statsmodels.api as sm # regression
from statsmodels.stats.outliers_influence import variance_inflation_factor # VIF test

# plt 한글깨짐 방지
font_name=font_manager.FontProperties(fname="c:/windows/fonts/malgun.ttf").get_name()
rc('font',family=font_name)


#### 데이터 불러오기

#일별, 국가별 확진자수 데이터 - 출처 : WHO
covid = pd.read_csv("c:/data/covid_project/WHO-COVID-19-global-data.csv")
covid.head()
covid.info()

#2020년 10월 06일을 기준으로 누적 확진자 수와 국가별 gdp, 면적, iq, 삶의 질, 인구수, 출생 및 사망률, 국가 성장률, moral machine의 law, 시민의식
data=pd.read_csv("c:/data/covid_project/covid_analysis_github.csv",encoding='euc-kr')
data.head()
data.info()


#### 전처리 및 파생변수 생성
### 1. 중국 데이터 제거
data=data[data.country_code!='CN'] #중국은 확진자를 제대로 집계하지 않고, 코로나 종식을 선언하다. 집계의 신뢰성을 믿을 수 없기 떄문에 제거

#팬데믹이 선언되었지만, 확진자가 200명 미만인 나라는 소국이거나 집계의 신뢰성을 믿을 수 없어 제거
(data.cumulative_cases<200).sum()/len(data) #44국, 18%
data=data[data.cumulative_cases>=200]
data.reset_index(drop=True,inplace=True) #인덱스 초기화

### 2. 파생변수 생성
data["infection_per_population"]=data["cumulative_cases"]/data["population"] # 인구수 대비 누적 확진자, ***주요 변수***
data["death_per_population"]=data["cumulative_deaths"]/data["population"] # 인구수 대비 누적 사망자
data["gdp_per_capita"]=data["gdp"]/data["population"] # 1인당 GDP
data["density"]=data["population"]/data["area"] # 인구밀도


### 3. 컬럼별 NA의 개수 확인
pd.DataFrame({'NA의 개수':pd.isnull(data).sum(),'전체 데이터 중 비율':np.round(pd.isnull(data).sum()/len(data)*100)})
msno.matrix(df=data.iloc[:, :], color=(0.1, 0.6, 0.8),figsize=(15,8),fontsize=7) # NA시각화

### 4. covid데이터 정제
covid.columns =covid.columns.str.lstrip().str.lower() # 컬럼이름 공백 제거, 소문자화
covid['date_reported']=pd.to_datetime(covid['date_reported']) # 날짜형으로 변환
covid.loc[pd.isnull(covid["country_code"]),"country_code"] ='NAM' #naibma의 약자가 NA라서 3약자인 NAM으로 수정


#### 분석 1.  zero-order correlation test
plt.figure(figsize=(16,13))
sns.heatmap(data=data.corr(),annot=True, fmt='.2f',linewidths=.5,cmap="Blues")
# data.corr().to_csv("c:/data/covid_project/result1.csv")


# 인구당 누적 확진자에 대해 n,r,p-value table 생성
corr_table_infection_per_population=pd.DataFrame({'n':[],'r':[],'p-val':[]})
for i in range(len(data.columns)):
    try :
        a=pg.corr(x=data['infection_per_population'], y=data[data.columns[i]])[["n","r","p-val"]]
        a.index = [data.columns[i]]
        a=np.round(a,3)
        corr_table_infection_per_population=pd.concat([corr_table_infection_per_population,a])
    except :
        pass # 문자형 변수는 pass
corr_table_infection_per_population['star']=(corr_table_infection_per_population['p-val']<=0.05).apply(lambda x:"*" if x else "")
corr_table_infection_per_population # 인구당 누적 확진자와 다른 변수들간의 correlation table


#### 분석 2. 코로나 초기 확산 속도 구하고 변수들간의 상관관계 파악
'''
국가별 코로나 첫 확진자가 나온 시기부터 10.06일 까지 경과된 일수를 표시한다.
그리고 나서 누적확진자의 10%가 되는 시점을 '코로나의 초기 확상 속도'라고 정의하고 분석 시행1
'''


## 1. 국가별로 임의의 일련번호 매기기
codes=pd.DataFrame({"country_code":data.country_code.unique(),"code":range(1,len(data.country_code.unique())+1)})
codes

covid=pd.merge(covid,codes) #일별 확진자 데이터와 join
covid.info()

data=pd.merge(data,codes) # 10월 06일누적 확진자 데이터와 join

## 2. 국가별 코로나가 처음 발생한 날짜로 부터 현재가지 경과일을 구하고 join

# 누적확진자가 0인 일수는 삭제
covid=covid[covid.cumulative_cases !=0]

# 첫환자 발생시점부터 현재까지 경과일 구하기
a=[]
for i in list(covid.code.unique()):
    list(range(1,len(covid[covid['code']==i])+1))
    a+=list(range(1,len(covid[covid['code']==i])+1))

covid["day"]=a # 일별 확진자 수 데이터에 추가

## 3. 코로나 확진자의 전날대비 확진자 수를 구하기
a=[]
for i in list(covid.code.unique()):
    b=covid.loc[covid.code == i, "cumulative_cases"].pct_change() * 100
    a+= list(b)

covid["increase"] = a

##각 날짜별로 총 누적 확진자 수 대비 당일 누적 확진자 계산
a=[]
for i in list(covid.code.unique()):
    for j in covid.loc[covid.code == i, "cumulative_cases"].tolist():
        b=j/covid.loc[covid.code==i,"cumulative_cases"].max()*100
        a.append(b)

covid["speed"]=a

covid.info()
covid.head(50)

## 총 누적 확진자수 대비 당일 누적확진자 수가 10%가 되는 시점 구하기
table_10=covid.loc[covid.speed<11,["speed","code","day"]].groupby("code").max()

table_10=table_10.reset_index()
table_10.columns = ["code","speed","10day"]

corona=pd.merge(data,table_10[["code","10day"]],left_on="code",right_on="code")


#시각화
plt.hist(corona["10day"],color="grey")
plt.xlabel("경과 일수",size=14)
plt.ylabel("Count",size=13)
plt.ylim(0,40)


### 분석 2-1. 10%가 되는 시점(초기 코로나 확산 시점)과 correlation
plt.figure(figsize=(15,15))
sns.heatmap(data=corona.corr(),annot=True, fmt='.2f',linewidths=.5,cmap="Blues")
# corona.corr().to_csv("c:/data/covid_project/result21.csv")

corr_table_10day=pd.DataFrame({'n':[],'r':[],'p-val':[]})
for i in range(len(corona.columns)):
    try :
        a=pg.corr(x=corona['10day'], y=corona[corona.columns[i]])[["n","r","p-val"]]
        a.index = [corona.columns[i]]
        a=np.round(a,3)
        corr_table_10day=pd.concat([corr_table_10day,a])
    except :
        pass
corr_table_10day['star']=(corr_table_10day['p-val']<=0.05).apply(lambda x:"*" if x else "")
corr_table_10day
# corona.to_csv("c:/data/corona_spss.csv",index=False)


#### 분석 2-2.regression test
feature=['10day','iq','quality_of_life_total', 'citizenship','gdp_per_capita']

regression_table =corona[feature]
regression_table=regression_table.dropna()

model=sm.OLS(regression_table['10day'],regression_table.iloc[:,1:10]).fit()
model.summary()

##다중 공선성 체크
regression_table[["quality_of_life_total","citizenship"]].corr() #높은 상관관계가 관측됨

# 모델의 VIF test
feature_independent=['iq','quality_of_life_total', 'citizenship','gdp_per_capita']
regression_table =corona[feature_independent]
regression_table=regression_table.dropna()
vif = pd.DataFrame()
vif["VIF Factor"] = [variance_inflation_factor(regression_table.values, i) for i in range(regression_table.shape[1])]
vif["features"] = regression_table.columns
vif

# VIF 계수가 가장 큰 'quality_of_life_total'을 제거하고 다시 VIF test
feature_independent=['iq', 'citizenship','gdp_per_capita'] # 삶의 질 점수를 제외하고 regression
regression_table =corona[feature_independent]
regression_table=regression_table.dropna()
vif = pd.DataFrame()
vif["VIF Factor"] = [variance_inflation_factor(regression_table.values, i) for i in range(regression_table.shape[1])]
vif["features"] = regression_table.columns
vif


# VIF 계수가 가장 큰 'citizenship'을 제거하고 다시 regression test
feature=['10day', 'iq','gdp_per_capita'] # 삶의 질 점수를 제외하고 regression
regression_table =corona[feature]
regression_table=regression_table.dropna()
model=sm.OLS(regression_table['10day'],regression_table.iloc[:,1:10]).fit()
model.summary()


#'10day'와 상관관계가 있던 일부 변수들과의 sactter plot
sns.set(style='white', font_scale=1.2) #세팅

#moral_machine_law와의 scatter plot
g = sns.JointGrid(data=corona,  x='moral_machine_law', y='10day')
g = g.plot_joint(sns.regplot, color="xkcd:muted blue")
g.set_axis_labels("obey the law","코로나 확산 속도")
g = g.plot_marginals(sns.distplot, kde=False, bins=12, color="xkcd:bluey grey")
plt.tight_layout()


# quality_of_life_health와의 scatter plot
g = sns.JointGrid(data=corona, x='10day', y='quality_of_life_health')
g = g.plot_joint(sns.regplot, color="xkcd:muted blue")
g.set_axis_labels("코로나 확산 속도", "qualityof life(health)")
g = g.plot_marginals(sns.distplot, kde=False, bins=12, color="xkcd:bluey grey")
plt.tight_layout()

# iq와의 scatter plot
g = sns.JointGrid(data=corona, x='10day', y='iq')
g = g.plot_joint(sns.regplot, color="xkcd:muted blue")
g.set_axis_labels("코로나 확산 속도", "IQ score")
g = g.plot_marginals(sns.distplot, kde=False, bins=12, color="xkcd:bluey grey")
plt.tight_layout()

# gdp_per_capita와의 scatter plot
g = sns.JointGrid(data=corona, x='10day', y='gdp_per_capita')
g = g.plot_joint(sns.regplot, color="xkcd:muted blue")
g = g.plot_marginals(sns.distplot, kde=False, bins=12, color="xkcd:bluey grey")
plt.tight_layout()



### 분석 2-3. 국가의 기온과 1인당 gdp의 상화작용,조절효과 검증

inter=corona[["code","10day","temperature","gdp_per_capita"]]
inter=inter.dropna()

inter['temperature_z']=ss.zscore(inter['temperature']) #z_score로 변경
inter['gdp_per_capita_z']=ss.zscore(inter['gdp_per_capita'])

inter["interaction"]=inter['temperature_z']*inter['gdp_per_capita_z'] #interaction term 생성

# interaction term을 넣은 모델과 넣지 않은 모델의 변화된 R-square 값을 확인
model1=sm.OLS(inter['10day'],inter.iloc[:,4:6]).fit()
model1.summary()

model2=sm.OLS(inter['10day'],inter.iloc[:,4:]).fit()
model2.summary() #interaction term의 p-value는 유의한 수준

model2.rsquared-model1.rsquared #13%의 변화량이 발생


## 조절효과의시각화

# 기온,1인당 gdp이 평균보다 높은면 1 낮으면 0로 코딩
inter['temperature_2']=inter['temperature_z'].apply(lambda x: 1 if x>=0 else 0)
inter['gdp_per_capita_2']=inter['gdp_per_capita_z'].apply(lambda x: 1 if x>=0 else 0)

pd.crosstab(inter['temperature_2'],inter['gdp_per_capita_2'])
pd.crosstab(inter['temperature_2'],inter['gdp_per_capita_2'], values=inter['10day'], aggfunc='mean')

# 추가적인 분석을 위해 2변수의 평균값을 확인
inter['gdp_per_capita'].mean() #18039
inter['temperature'].mean() # 23


# 두 변수의 평균과 interaction effect의 분석결과를 기준으로 4개의 집단으로 나눠 그래프 그리기
inter['group']=inter['temperature_2'].map(str)+inter['gdp_per_capita_2'].map(str)
inter['group'].value_counts()

plt_table=inter[["group","10day"]].groupby("group").mean()
plt_table["temperature"]=["low","low",'high','high']
plt_table["gdp_per_capita"]=["low",'high',"low",'high']

plt.plot(plt_table.loc[plt_table['gdp_per_capita']=="low",'temperature'],plt_table.loc[plt_table['gdp_per_capita']=="low",'10day'],label ="low_gdp_per_capita",color="r",linestyle="--")
plt.plot(plt_table.loc[plt_table['gdp_per_capita']=="high",'temperature'],plt_table.loc[plt_table['gdp_per_capita']=="high",'10day'],label ="high_gdp_per_capita",color="b",linestyle=":")
plt.legend()
plt.xlabel("temperature")
plt.ylabel("코로나 확산 속도")
plt.ylim(40,100)

# gdp_per_capita와 기온의의 scatter plot
g = sns.JointGrid(data=corona, x='temperature', y='gdp_per_capita')
g = g.plot_joint(sns.regplot, color="xkcd:muted blue")
g.set_axis_labels("기온", "1인당GDP")
g = g.plot_marginals(sns.distplot, kde=False, bins=12, color="xkcd:bluey grey")
plt.tight_layout()

