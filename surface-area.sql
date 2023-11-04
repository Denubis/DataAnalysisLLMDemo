sqlite> select aswkt(st_expand(MakePoint(20.5, -64.5, 4326),0.5));
+----------------------------------------------------+
| aswkt(st_expand(MakePoint(20.5, -64.5, 4326),0.5)) |
+----------------------------------------------------+
| POLYGON((20 -65,21 -65,21 -64,20 -64,20 -65))      |
+----------------------------------------------------+
sqlite> select area(st_expand(MakePoint(20.5, -64.5, 4326),0.5),true);
+--------------------------------------------------------+
| area(st_expand(MakePoint(20.5, -64.5, 4326),0.5),true) |
+--------------------------------------------------------+
| 5357200336.53382                                       |
+--------------------------------------------------------+
create table surface_area(longitude real, latitude real, area real);

insert into surface_area (longitude, latitude, area) SELECT distinct longitude, latitude, area(st_expand(MakePoint(longitude, latitude, 4326),0.5),true) from argo_data;
