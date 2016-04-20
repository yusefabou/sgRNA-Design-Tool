from django.conf.urls import include, url
from django.contrib import admin

urlpatterns = [
    url(r'^crispr_app/', include('crispr_app.urls')),
    url(r'^admin/', include(admin.site.urls)),
]
