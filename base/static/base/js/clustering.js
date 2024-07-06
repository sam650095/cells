let methodtextnode;
document.addEventListener("DOMContentLoaded", async function () {
    const csrftoken = getCookie("csrftoken");
    const preload_clustering_results = await fetchAPI("/api/preloadclustering", 0, csrftoken);
    console.log(preload_clustering_results);
    show_clustering_method(preload_clustering_results);
  });
// show selection
function show_clustering_method(preload_clustering_results) {
    const clustering_method = document.getElementById('clustering_method');
    const clustering_method_ul = document.getElementById('preloadul');
    
    clustering_method_ul.innerHTML = '';
    methodtextnode = document.createTextNode('none');
    clustering_method.insertBefore(methodtextnode, clustering_method.firstChild);
    for (let i = 0; i < preload_clustering_results.methods.length; i++) {
        const li = document.createElement('li');
        li.innerHTML = `
                <a onclick="methodtextnodefix('${preload_clustering_results.methods[i]}')"
                class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2"
            >${preload_clustering_results.methods[i]}</a>
        `;
        clustering_method_ul.appendChild(li);
    }
}
function methodtextnodefix(newText){
    methodtextnode.nodeValue = newText;
  }
// fetch api
async function fetchAPI(url, formData, csrftoken) {
    const response = await fetch(url, {
      method: "POST",
      body: formData,
      headers: {
        "X-CSRFToken": csrftoken,
      },
    });
  
    if (!response.ok) {
      throw new Error("Network response was not ok");
    }
  
    return response.json();
  }
  