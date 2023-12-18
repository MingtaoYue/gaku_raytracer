# ビルド
```
git clone https://github.com/MingtaoYue/gaku_raytracer.git
cd gaku_raytracer
g++ -o raytracer raytracer.cpp
./raytracer 2
```
数字のパラメータはサブピクセル毎のサンプリング光線数です、デフォルトは1で、大きくなると、よりレンダリング質が良くなりますが、時間がかかる。
# 説明
パストレーシング手法を用いて実装したレイトレーシングソフトウェアです。
