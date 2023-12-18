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
- ピクセル毎に複数の光線を追跡し、平均ラディアンスを算出してレンダリング
- レンダリング方程式ではモンテカルロ法を用いたサンプリングにより近似解を得る
- 光線ごとに、シーン内の物体との交差点を再帰的に計算し、反射または屈折を模擬
- ロシアンルーレット法により、計算の再帰を確率的に打ち切り、真値に対する期待値と一致する結果を効率的に近似させる
# レンダリング例
サンプリング光線数=1
<img width="900" alt="Screenshot 2023-12-18 at 11 31 24" src="https://github.com/MingtaoYue/gaku_raytracer/assets/127390549/f77e6170-b072-4eda-8bc5-7ecd109e70a0">

サンプリング光線数=200
<img width="900" alt="Screenshot 2023-12-18 at 11 33 14" src="https://github.com/MingtaoYue/gaku_raytracer/assets/127390549/a4146cbb-5048-45f1-b060-89d619ab3c78">



